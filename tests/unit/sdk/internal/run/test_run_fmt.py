"""
@author: cunyue
@file: test_run_helper.py
@time: 2026/3/12 16:54
@description: 测试run工具函数
"""

import pytest

import swanlab.sdk.internal.run.fmt as fmt

# ─────────────────────────── flatten_dict ────────────────────────────


class TestFlattenDict:
    """测试 fmt.flatten_dict 函数"""

    def test_empty_dict(self):
        """空字典直接返回空字典"""
        assert fmt.flatten_dict({}) == {}

    def test_flat_dict_unchanged(self):
        """无嵌套字典不做任何变换"""
        assert fmt.flatten_dict({"a": 1, "b": "x"}) == {"a": 1, "b": "x"}

    def test_nested_one_level(self):
        """单层嵌套展开"""
        assert fmt.flatten_dict({"a": {"b": 1}}) == {"a/b": 1}

    def test_nested_deep(self):
        """深层嵌套展开"""
        assert fmt.flatten_dict({"a": {"b": {"c": 1}}}) == {"a/b/c": 1}

    def test_mixed_flat_and_nested(self):
        """平铺与嵌套混合"""
        result = fmt.flatten_dict({"x": 0, "a": {"b": 1, "c": 2}})
        assert result == {"x": 0, "a/b": 1, "a/c": 2}

    def test_int_key_converted_to_str(self):
        """整数 key 自动转为字符串"""
        assert fmt.flatten_dict({1: {2: "v"}}) == {"1/2": "v"}  # type: ignore

    def test_duplicate_key_latter_wins(self, monkeypatch):
        """出现同名展开键时，后者覆盖前者，并触发一次警告"""
        mock_warning = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: mock_warning.append(args))

        data = {"a": {"b": 1}, "a/b": 99}
        result = fmt.flatten_dict(data)

        assert len(mock_warning) == 1
        assert result["a/b"] == 99

    def test_returns_parent_dict_when_provided(self):
        """传入 parent_dict 时结果合并到其中"""
        parent = {"existing": "value"}
        result = fmt.flatten_dict({"a": 1}, parent_dict=parent)
        assert result is parent
        assert result == {"existing": "value", "a": 1}


# ─────────────────────────── validate_key ────────────────────────────


class TestValidateKey:
    """测试 fmt.validate_key 函数"""

    @pytest.fixture(autouse=True)
    def reset_warned_keys(self):
        """每个测试前清空全局警告缓存，避免用例间相互污染"""
        fmt._WARNED_KEYS.clear()
        yield
        fmt._WARNED_KEYS.clear()

    @pytest.mark.parametrize("key", ["loss", "train/acc", "v1.0_metric", "loss * 1000", "my key", "a@b#c"])
    def test_valid_key_unchanged(self, key):
        """合法 key 原样返回（内部字符宽松接受）"""
        assert fmt.validate_key(key) == key

    def test_strip_leading_trailing(self, monkeypatch):
        """头尾的空白、. 和 / 被剥离"""
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: None)
        assert fmt.validate_key("  loss  ") == "loss"
        assert fmt.validate_key("./key/") == "key"
        assert fmt.validate_key("...key...") == "key"

    @pytest.mark.parametrize(
        "invalid, expected",
        [
            ("my key", "my key"),  # 内部空格被接受，原样返回
            ("a@b#c", "a@b#c"),  # 内部特殊字符被接受，原样返回
        ],
    )
    def test_internal_chars_accepted(self, invalid, expected):
        """内部字符只要不是控制字符，都原样保留"""
        assert fmt.validate_key(invalid) == expected

    def test_unicode_word_chars_allowed(self):
        """Python \\w 匹配 Unicode 字符（汉字等），不会被替换"""
        # \w 在 Python3 默认模式下匹配 Unicode 字母数字，汉字合法
        result = fmt.validate_key("中文key")
        assert result == "中文key"

    def test_truncate_long_key(self):
        """超长 key 被截断到 max_len"""
        long_key = "a" * 300
        result = fmt.validate_key(long_key, max_len=255)
        assert len(result) == 255
        assert result == "a" * 255

    def test_custom_max_len(self, monkeypatch):
        """自定义 max_len 生效"""
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: None)
        assert fmt.validate_key("abcdef", max_len=3) == "abc"

    @pytest.mark.parametrize("key, expected", [(42, "42"), (3.14, "3.14")])
    def test_non_string_input(self, key, expected):
        """非字符串 key 被转为字符串处理"""
        assert fmt.validate_key(key) == expected  # type: ignore

    def test_empty_after_sanitization_raises(self):
        """清洗后为空时抛出 ValueError"""
        with pytest.raises(ValueError, match="invalid or empty after sanitization"):
            fmt.validate_key("   ")
        with pytest.raises(ValueError, match="invalid or empty after sanitization"):
            fmt.validate_key("...")

    def test_warning_fires_once_per_key(self, monkeypatch):
        """相同 key 被头尾剥离时，只触发一次警告"""
        mock_warning = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: mock_warning.append(args))

        fmt.validate_key("  bad_key  ")
        fmt.validate_key("  bad_key  ")
        fmt.validate_key("  bad_key  ")

        assert len(mock_warning) == 1

    def test_warning_fires_for_different_keys(self, monkeypatch):
        """不同被剥离的 key 各自独立触发一次警告"""
        mock_warning = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: mock_warning.append(args))

        fmt.validate_key("  bad_key1  ")
        fmt.validate_key("  bad_key2  ")

        assert len(mock_warning) == 2


# ─────────────────────────── save validators ────────────────────────────


class TestSaveValidators:
    """测试 save 相关校验函数"""

    @pytest.mark.parametrize("policy", ["now", "end", "live"])
    def test_safe_validate_save_policy_accepts_valid_policy(self, policy):
        assert fmt.safe_validate_save_policy(policy) == policy

    @pytest.mark.parametrize(("policy", "expected"), [("NOW", "now"), ("End", "end"), ("LIVE", "live")])
    def test_safe_validate_save_policy_normalizes_case(self, policy, expected):
        assert fmt.safe_validate_save_policy(policy) == expected

    @pytest.mark.parametrize("policy", ["invalid", "", None, 1])
    def test_safe_validate_save_policy_rejects_invalid_policy(self, policy):
        assert fmt.safe_validate_save_policy(policy) is None


class TestResolveSavePaths:
    """测试 save 路径解析与校验"""

    def test_resolve_relative_glob_uses_cwd_as_base(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        resolved = fmt.resolve_save_paths("model.pt")

        assert resolved == ((tmp_path / "model.pt").resolve(), tmp_path.resolve())

    @pytest.mark.parametrize("glob_str", [None, 1, object()])
    def test_invalid_glob_type_returns_none(self, glob_str, monkeypatch):
        errors = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.error", lambda *args: errors.append(args))

        assert fmt.resolve_save_paths(glob_str) is None  # type: ignore[arg-type]
        assert errors

    def test_invalid_glob_bytes_returns_none(self, monkeypatch):
        errors = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.error", lambda *args: errors.append(args))

        assert fmt.resolve_save_paths(b"\xff") is None
        assert errors

    @pytest.mark.parametrize("base_path", [1, object(), []])
    def test_invalid_base_path_type_returns_none(self, base_path, monkeypatch):
        errors = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.error", lambda *args: errors.append(args))

        assert fmt.resolve_save_paths("model.pt", base_path=base_path) is None  # type: ignore[arg-type]
        assert errors

    def test_glob_outside_base_path_returns_none(self, tmp_path, monkeypatch):
        errors = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.error", lambda *args: errors.append(args))
        base_path = tmp_path / "base"
        base_path.mkdir()

        assert fmt.resolve_save_paths(str(tmp_path / "model.pt"), base_path=base_path) is None
        assert errors

    @pytest.mark.parametrize("glob_str", ["s3://bucket/model.pt", "gs://bucket/model.pt"])
    def test_cloud_storage_url_returns_none(self, glob_str, monkeypatch):
        warnings = []
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.warning", lambda *args: warnings.append(args))

        assert fmt.resolve_save_paths(glob_str) is None
        assert warnings

    def test_path_resolve_error_returns_none(self, monkeypatch):
        errors = []
        original_resolve = fmt.Path.resolve

        def raise_on_model_path(path):
            if str(path) == "model.pt":
                raise OSError("failed to resolve")
            return original_resolve(path)

        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.console.error", lambda *args: errors.append(args))
        monkeypatch.setattr("swanlab.sdk.internal.run.fmt.Path.resolve", raise_on_model_path)

        assert fmt.resolve_save_paths("model.pt") is None
        assert errors
