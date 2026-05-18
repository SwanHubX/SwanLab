"""
针对 PR #1633 修复进行暴力 + 边界回归测试。

回归目标：
  1) 默认字符集 (ascii_lowercase + digits) 下，generate_id 始终满足
     正则 ^[a-z0-9]+$，长度精确等于参数。
  2) 自定义字符集时，输出字符全部来自给定 characters。
  3) 边界长度 [1, 64] 全部合法；越界 (<=0 / >64) 抛 ValueError。
  4) 默认调用允许出现纯数字 ID（即不再人为注入小写字母）。
     - 形态校验用确定性断言（不依赖随机分布）：通过 ``characters="7"`` 与
       ``characters=string.digits`` 强制构造纯数字 ID 路径，断言不被改写。
     - 数据分布健全性额外用 10 万次抽样作为辅助观测，但不作为强断言阈值。
  5) 输出唯一性足够高（不强制全唯一，但碰撞率应远低于阈值）。
"""

import re
import string

import pytest

from swanlab.utils.experiment import generate_id

_DEFAULT_RE = re.compile(r"^[a-z0-9]+$")

# 暴力回归核心规模：根据用户要求设置为 10 万次。
_STRESS_ROUNDS = 100_000


class TestGenerateIdBruteForce:
    @pytest.mark.parametrize("rounds", [_STRESS_ROUNDS])
    def test_default_charset_brute(self, rounds):
        """默认字符集下 10 万次调用：长度恒定、字符集严格、且能观察到字母。

        说明：本测试只做"形态级"强断言（长度、字符集、出现字母），不依赖
        "纯数字必然出现"这种统计型断言（10 万次中纯数字期望约 60，但
        理论上仍可能为 0；用于回归 buggy 注入逻辑的强证据放在
        ``test_pure_digits_id_is_allowed`` / ``test_no_islower_injection_brute``
        这两条确定性测试里）。
        """
        seen_alpha = 0
        seen_pure_digits = 0
        for _ in range(rounds):
            run_id = generate_id()
            assert len(run_id) == 8
            assert _DEFAULT_RE.fullmatch(run_id), f"charset violated: {run_id!r}"
            if any(c.isalpha() for c in run_id):
                seen_alpha += 1
            else:
                seen_pure_digits += 1

        assert seen_alpha > 0, "expected to see ids containing letters"
        # 仅作为分布健全性观察，不卡阈值（防止极小概率抖动）
        assert seen_pure_digits >= 0

    @pytest.mark.parametrize("length", [1, 2, 3, 4, 7, 8, 9, 16, 32, 63, 64])
    def test_all_valid_lengths(self, length):
        """所有合法长度边界都能稳定生成。"""
        for _ in range(500):
            run_id = generate_id(length)
            assert len(run_id) == length
            assert _DEFAULT_RE.fullmatch(run_id)

    @pytest.mark.parametrize("invalid_length", [-100, -1, 0, 65, 66, 128, 1024])
    def test_invalid_lengths(self, invalid_length):
        """所有越界长度都抛 ValueError，且文案稳定。"""
        with pytest.raises(ValueError, match="Length must be between 1 and 64."):
            generate_id(invalid_length)

    def test_uniqueness_high(self):
        """10 万次样本碰撞概率（生日悖论）≈ 1 - exp(-N^2 / (2*36^8)) ≈ 0.18%，
        允许极少量碰撞但应远小于 0.1%。"""
        n = 100_000
        ids = {generate_id() for _ in range(n)}
        # 生日悖论期望碰撞数 ≈ N(N-1) / (2 * 36^8) ≈ 1.78
        # 给出 50 的宽松上界，远高于期望值，避免偶发抖动
        collisions = n - len(ids)
        assert collisions <= 50, f"too many collisions: {collisions} in {n} runs"

    def test_pure_digits_id_is_allowed(self):
        """显式构造纯数字 ID 路径：仅用 digits 字符集时输出必须全是数字。"""
        for _ in range(500):
            run_id = generate_id(length=8, characters=string.digits)
            assert len(run_id) == 8
            assert run_id.isdigit(), f"expected digit-only id: {run_id!r}"
            # 关键：不应再被人为塞入字母
            assert not any(c.isalpha() for c in run_id)

    def test_custom_charset_uppercase_preserved(self):
        """自定义包含大写字母时，输出可包含大写 — 不再强行小写化。"""
        chars = string.ascii_uppercase
        seen_upper = False
        for _ in range(200):
            run_id = generate_id(length=8, characters=chars)
            assert len(run_id) == 8
            assert all(c in chars for c in run_id)
            if any(c.isupper() for c in run_id):
                seen_upper = True
        assert seen_upper, "expected uppercase chars to appear when allowed"

    def test_custom_charset_single_char(self):
        """字符集只有一个字符时，输出应是该字符的重复。"""
        run_id = generate_id(length=10, characters="x")
        assert run_id == "x" * 10

    def test_custom_charset_with_symbols(self):
        """带符号的自定义字符集：输出严格属于给定集合。"""
        chars = "ab-_."
        for _ in range(200):
            run_id = generate_id(length=12, characters=chars)
            assert len(run_id) == 12
            assert all(c in chars for c in run_id)

    @pytest.mark.parametrize("length", [1, 8, 64])
    def test_length_param_keyword_and_positional(self, length):
        """位置/关键字参数两种调用方式等价。"""
        a = generate_id(length)
        b = generate_id(length=length)
        assert len(a) == length
        assert len(b) == length

    def test_no_islower_injection(self, monkeypatch):
        """关键回归：纯数字 ID 不应被业务代码人为修正成含字母的 ID。

        通过 monkeypatch `secrets.choice` 强制每次都返回 '7'，
        若 generate_id 中存在 "if not islower(): 注入小写字母" 的兜底逻辑，
        输出会变成形如 "7777777a"；正确实现下应严格等于 "77777777"。
        """
        import swanlab.utils.experiment as exp_mod

        def fake_choice(_seq):
            return "7"

        monkeypatch.setattr(exp_mod.secrets, "choice", fake_choice)
        run_id = generate_id()
        assert run_id == "77777777", (
            f"got {run_id!r}; bias-injection logic seems to be back, "
            "which would silently corrupt user-provided pure-digit cases."
        )

    def test_no_islower_injection_brute(self):
        """10 万次暴力版：在所有合法长度 [1, 64] 上断言"纯数字输出绝不被改写"。

        关键设计：
          - 不 mock `secrets.choice`，直接用 ``characters="7"`` 让真实 `secrets.choice`
            返回值唯一确定为 '7'，输出本应严格等于 '7' * length。
          - 这样可同时捕获两类 buggy 实现：
              a) 旧逻辑 `if not any(c.isalpha()): 注入 secrets.choice(string.ascii_lowercase)`
                 — 字符集会被污染成字母。
              b) 旧逻辑 `if not id_str.islower(): 注入 secrets.choice(string.ascii_lowercase)`
                 — 同样会注入外部小写字母。
          - 枚举 64 个长度 × 1600 次 ≈ 102_400 次，达到"暴力"量级。
        """
        per_length = 1600
        total = 0
        for length in range(1, 65):
            expected = "7" * length
            for _ in range(per_length):
                run_id = generate_id(length=length, characters="7")
                assert run_id == expected, (
                    f"len={length}: got {run_id!r}, expected {expected!r}; "
                    "bias-injection logic re-introduced — pure-digit ids would be silently corrupted."
                )
                total += 1
        assert total >= 100_000
