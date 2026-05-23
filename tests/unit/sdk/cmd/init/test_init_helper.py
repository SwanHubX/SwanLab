"""
@author: cunyue
@file: test_init_helper.py
@time: 2026/3/10 20:52
@description: 测试 init 中的辅助函数
"""

import json
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.cmd.init import (
    compatible_kwargs,
    ensure_run_dir,
    load_config,
    prompt_init_mode,
    set_nested_value,
)


# ==========================================
# Fixtures
# ==========================================
@pytest.fixture
def mock_settings():
    """构造一个轻量级的 Settings Mock 对象"""
    settings = MagicMock()
    settings.run.config = None
    settings.mode = "online"
    settings.interactive = True
    settings.web_host = "https://swanlab.cn"
    return settings


# ==========================================
# Tests for `load_config`
# ==========================================
def test_load_config_basic(mock_settings):
    """测试基础类型的加载 (None 和 dict)"""
    assert load_config(mock_settings, None) == {}
    assert load_config(mock_settings, {"key": "value"}) == {"key": "value"}


def test_load_config_files(mock_settings, tmp_path):
    """测试从 JSON 和 YAML 文件中正确加载配置"""
    # 测试 JSON
    json_file = tmp_path / "config.json"
    json_file.write_text(json.dumps({"a": 1}))
    assert load_config(mock_settings, str(json_file)) == {"a": 1}

    # 测试 YAML
    yaml_file = tmp_path / "config.yaml"
    yaml_file.write_text("b: 2\nc: 3")
    assert load_config(mock_settings, yaml_file) == {"b": 2, "c": 3}


def test_load_config_exceptions(mock_settings):
    """测试异常分支：文件不存在和类型错误"""
    with pytest.raises(FileNotFoundError):
        load_config(mock_settings, "not_exist_file.json")

    with pytest.raises(ValueError, match="Invalid config type"):
        load_config(mock_settings, 123)  # type: ignore


# ==========================================
# Tests for `prompt_init_mode`
# ==========================================
def test_prompt_fast_exits(mock_settings, monkeypatch):
    """测试快速跳过的条件：非 online 模式、非交互模式、已登录"""
    # 1. 模拟已登录 (client.exists 返回 True)
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: True)
    assert prompt_init_mode(mock_settings) == "online"

    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)

    # 2. 模拟非交互模式
    mock_settings.interactive = False
    assert prompt_init_mode(mock_settings) == "online"

    # 3. 模拟非云端模式
    mock_settings.interactive = True
    mock_settings.mode = "local"
    assert prompt_init_mode(mock_settings) == "local"


def test_prompt_auto_login(mock_settings, monkeypatch):
    """测试本地已存在 apikey 时直接返回 online 模式"""
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)
    # 设置 api_key 模拟已存在凭证
    mock_settings.api_key = "fake-key"

    mock_login_raw = MagicMock(return_value=True)
    monkeypatch.setattr("swanlab.sdk.cmd.init.login_raw", mock_login_raw)

    assert prompt_init_mode(mock_settings) == "online"
    # 有 api_key 时不再调用 login_raw，登录推迟到后续流程


@pytest.mark.parametrize(
    "inputs, mock_login_success, expected_mode",
    [
        # 直接选 3 -> 切换离线模式
        (["3"], False, "offline"),
        # 直接选 1 -> 触发交互登录并返回成功
        (["1"], True, "online"),
        # 直接选 2 -> 触发交互登录并返回失败
        (["2"], False, "online"),
        # 乱输一通后选 3 -> 循环容错测试
        (["invalid", "wrong", "3"], False, "offline"),
    ],
)
def test_prompt_interactive_choices(mock_settings, monkeypatch, inputs, mock_login_success, expected_mode):
    """使用参数化和猴子补丁极致压缩终端交互的测试代码"""
    monkeypatch.setattr("swanlab.sdk.cmd.init.client.exists", lambda: False)
    mock_settings.api_key = None

    # 模拟 login_raw 返回值
    monkeypatch.setattr("swanlab.sdk.cmd.init.login_cli", lambda **kwargs: mock_login_success)

    # 模拟 input，通过迭代器按顺序弹出输入值
    input_iterator = iter(inputs)
    monkeypatch.setattr("builtins.input", lambda prompt: next(input_iterator))

    assert prompt_init_mode(mock_settings) == expected_mode


class TestSetNestedValue:
    """测试 set_nested_value 函数"""

    def test_simple_key(self):
        """测试简单的单层键"""
        d = {}
        set_nested_value(d, "a", 1)
        assert d == {"a": 1}

    def test_nested_key(self):
        """测试嵌套键路径"""
        d = {}
        set_nested_value(d, "a.b.c", 123)
        assert d == {"a": {"b": {"c": 123}}}

    def test_update_existing_nested_value(self):
        """测试更新已存在的嵌套值"""
        d = {"a": {"b": {"c": 1}}}
        set_nested_value(d, "a.b.c", 999)
        assert d == {"a": {"b": {"c": 999}}}

    def test_partial_existing_path(self):
        """测试部分路径已存在的情况"""
        d = {"a": {"x": 1}}
        set_nested_value(d, "a.b.c", 2)
        assert d == {"a": {"x": 1, "b": {"c": 2}}}

    def test_none_value_not_set(self):
        """测试 None 值不会被设置"""
        d = {}
        set_nested_value(d, "a", None)
        assert d == {}

    def test_none_value_nested_not_set(self):
        """测试嵌套路径的 None 值不会被设置"""
        d = {"x": 1}
        set_nested_value(d, "a.b.c", None)
        assert d == {"x": 1}


class TestCompatibleKwargs:
    """测试 compatible_kwargs 的参数兼容性映射"""

    def test_experiment_name_maps_to_nested_key(self):
        """experiment_name → experiment.name"""
        result = compatible_kwargs({}, experiment_name="my_exp")

        assert result == {"experiment": {"name": "my_exp"}}

    def test_notes_maps_to_description(self):
        """notes → experiment.description"""
        result = compatible_kwargs({}, notes="my description")

        assert result == {"experiment": {"description": "my description"}}

    def test_both_kwargs_mapped_simultaneously(self):
        """experiment_name 和 notes 同时传入时，均应被正确映射"""
        result = compatible_kwargs({}, experiment_name="exp", notes="desc")

        assert result["experiment"]["name"] == "exp"
        assert result["experiment"]["description"] == "desc"

    def test_no_extra_kwargs_returns_model_dict_unchanged(self):
        """无额外 kwargs 时，model_dict 原样返回"""
        d = {"key": "value"}
        result = compatible_kwargs(d)

        assert result == {"key": "value"}

    def test_none_values_not_written_to_dict(self):
        """experiment_name=None 和 notes=None 不应写入任何键"""
        result = compatible_kwargs({}, experiment_name=None, notes=None)

        assert result == {}


class TestEnsureRunDir:
    """测试 ensure_run_dir 的原子化目录创建与冲突重试"""

    def test_creates_run_dir_when_no_conflict(self, tmp_path):
        """无冲突时，应成功创建 run_dir 并返回其路径"""
        run_dir = ensure_run_dir(tmp_path, "abc123", retry_interval=0.01)

        assert run_dir.exists()
        assert run_dir.is_dir()
        assert "abc123" in run_dir.name

    def test_retries_on_conflict(self, tmp_path, monkeypatch):
        """run_dir 已被占用时，应等待后用新时间戳重试直到成功"""
        call_count = 0
        original_generate = __import__(
            "swanlab.sdk.cmd.init", fromlist=["_generate_run_dir_name"]
        )._generate_run_dir_name
        this_run_dir = ""

        def mock_generate(run_id, max_length, parallel="none", hostname_max_len=64):
            nonlocal call_count, this_run_dir
            call_count += 1
            if call_count == 1:
                return "run-20260101_000000-conflict", False
            this_run_dir, truncated = original_generate(run_id, max_length, parallel, hostname_max_len)
            return this_run_dir, truncated

        monkeypatch.setattr("swanlab.sdk.cmd.init._generate_run_dir_name", mock_generate)

        # 预先创建冲突目录
        (tmp_path / "run-20260101_000000-conflict").mkdir()

        run_dir = ensure_run_dir(tmp_path, "abc123", retry_interval=0.01)

        assert run_dir.exists()
        assert run_dir.name == this_run_dir
        assert call_count >= 2

    def test_returns_path_under_log_dir(self, tmp_path):
        """返回的路径应在 log_dir 下"""
        run_dir = ensure_run_dir(tmp_path, "test-id", retry_interval=0.01)

        assert run_dir.parent == tmp_path

    def test_run_dir_name_format(self, tmp_path):
        """run_dir 名称应匹配 run-{timestamp}-{run_id} 格式"""
        run_dir = ensure_run_dir(tmp_path, "myrun42", retry_interval=0.01)

        name = run_dir.name
        assert name.startswith("run-")
        assert "myrun42" in name
        parts = name.split("-", 2)  # ["run", "20260418_123456", "myrun42"]
        assert len(parts) == 3
        timestamp = parts[1]
        assert len(timestamp) == 15  # YYYYMMDD_HHMMSS
        assert "_" in timestamp

    def test_raises_when_max_retries_exceeded(self, tmp_path, monkeypatch):
        """超过最大重试次数仍无法创建唯一目录时，应抛出 RuntimeError"""
        monkeypatch.setattr(
            "swanlab.sdk.cmd.init._generate_run_dir_name",
            lambda _run_id, _max_length, _parallel="none", _hostname_max_len=64: (
                "run-20260101_000000-conflict",
                False,
            ),
        )

        # 预先创建冲突目录
        (tmp_path / "run-20260101_000000-conflict").mkdir()

        with pytest.raises(RuntimeError, match="Failed to create a unique run directory after 2 attempts"):
            ensure_run_dir(tmp_path, "abc123", max_retries=2, retry_interval=0.01)

    def test_keeps_short_name(self, tmp_path):
        """未超过长度限制时，应保持原始格式"""
        run_dir = ensure_run_dir(tmp_path, "myrun42", retry_interval=0.01)

        name = run_dir.name
        assert name.startswith("run-")
        assert "myrun42" in name
        assert "..." not in name

    def test_truncates_long_name(self, tmp_path, monkeypatch):
        """超过长度限制时，应截断并在目录名中体现 ... 省略"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", MagicMock())

        run_dir = ensure_run_dir(tmp_path, "a" * 512, retry_interval=0.01, dir_max_length=80)

        assert len(run_dir.name.encode("utf-8")) <= 80
        assert "..." in run_dir.name

    def test_truncation_same_prefix(self, tmp_path, monkeypatch):
        """不同 run_id 截断后前缀相同时，... 省略格式可能相同"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", MagicMock())

        run_dir1 = ensure_run_dir(tmp_path, "a" * 511 + "1", retry_interval=0.01, dir_max_length=80)
        # 用不同的 log_dir 避免冲突
        log_dir2 = tmp_path / "second"
        log_dir2.mkdir()
        run_dir2 = ensure_run_dir(log_dir2, "a" * 511 + "2", retry_interval=0.01, dir_max_length=80)

        assert "..." in run_dir1.name
        assert "..." in run_dir2.name
        assert len(run_dir1.name.encode("utf-8")) <= 80
        assert len(run_dir2.name.encode("utf-8")) <= 80

    def test_warns_only_for_created_truncated_name(self, tmp_path, monkeypatch):
        """重试创建目录时，只对最终创建成功的截断目录告警"""
        warning = MagicMock()
        call_count = 0

        def mock_generate(_run_id, _max_length, _parallel="none", _hostname_max_len=64):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                return "run-20260101_000000-conflict", True
            return "run-20260101_000001-success", True

        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", warning)
        monkeypatch.setattr("swanlab.sdk.cmd.init._generate_run_dir_name", mock_generate)
        (tmp_path / "run-20260101_000000-conflict").mkdir()

        run_dir = ensure_run_dir(tmp_path, "a" * 512, max_retries=2, retry_interval=0.01, dir_max_length=80)

        assert run_dir.name == "run-20260101_000001-success"
        warning.assert_called_once()
        assert "run-20260101_000001-success" in warning.call_args.args[0]

    def test_does_not_warn_when_truncated_name_never_created(self, tmp_path, monkeypatch):
        """全部重试失败时，不应对未使用的截断目录告警"""
        warning = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", warning)
        monkeypatch.setattr(
            "swanlab.sdk.cmd.init._generate_run_dir_name",
            lambda _run_id, _max_length, _parallel="none", _hostname_max_len=64: (
                "run-20260101_000000-conflict",
                True,
            ),
        )
        (tmp_path / "run-20260101_000000-conflict").mkdir()

        with pytest.raises(RuntimeError, match="Failed to create a unique run directory after 2 attempts"):
            ensure_run_dir(tmp_path, "a" * 512, max_retries=2, retry_interval=0.01, dir_max_length=80)

        warning.assert_not_called()

    def test_with_dir_name_creates_once(self, tmp_path):
        """指定 dir_name 时直接使用，仅创建一次"""
        run_dir = ensure_run_dir(tmp_path, "abc123", dir_name="my-custom-dir")

        assert run_dir.exists()
        assert run_dir.name == "my-custom-dir"

    def test_with_dir_name_conflict_raises(self, tmp_path):
        """指定 dir_name 时目录已存在，直接抛错"""
        (tmp_path / "my-custom-dir").mkdir()

        with pytest.raises(FileExistsError):
            ensure_run_dir(tmp_path, "abc123", dir_name="my-custom-dir")


class TestEnsureRunDirShared:
    """测试 ensure_run_dir 在 shared 模式下的行为"""

    def test_shared_name_contains_hostname_and_pid(self, tmp_path, monkeypatch):
        """shared 模式目录名应包含 hostname 和 pid"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "myhost")
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 99999)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01)

        name = run_dir.name
        assert name.startswith("run-")
        assert "myhost" in name
        assert "99999" in name
        assert "run42" in name

    def test_shared_hostname_sanitize(self, tmp_path, monkeypatch):
        """hostname 中的非法字符应被替换为 _"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "node:01.cluster/local")
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01)

        name = run_dir.name
        assert ":" not in name
        assert "/" not in name
        assert "node_01_cluster_local" in name

    def test_shared_hostname_preserves_chinese(self, tmp_path, monkeypatch):
        """中文 hostname 应被保留"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "我的电脑")
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01)

        name = run_dir.name
        assert "我的电脑" in name

    def test_shared_hostname_all_illegal_fallback(self, tmp_path, monkeypatch):
        """全非法字符 hostname 应回退为 unknown_host"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "<>:|?*#%")
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01)

        name = run_dir.name
        assert "unknown_host" in name

    def test_shared_hostname_truncation(self, tmp_path, monkeypatch):
        """超长 hostname 应被截断为 abc...yz 格式"""
        long_hostname = "a" * 100
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: long_hostname)
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01)

        name = run_dir.name
        assert "..." in name
        assert len(name.encode("utf-8")) <= 255

    def test_shared_truncates_long_run_id(self, tmp_path, monkeypatch):
        """shared 模式下长 run_id 截断后不超限"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "myhost")
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", MagicMock())

        run_dir = ensure_run_dir(tmp_path, "a" * 512, parallel="shared", retry_interval=0.01, dir_max_length=80)

        assert len(run_dir.name.encode("utf-8")) <= 80
        assert "myhost" in run_dir.name
        assert "12345" in run_dir.name

    def test_none_mode_no_hostname(self, tmp_path):
        """parallel="none" 时目录名不应包含 hostname"""
        run_dir = ensure_run_dir(tmp_path, "run42", retry_interval=0.01)

        name = run_dir.name
        assert name.startswith("run-")
        parts = name.split("-")
        # none 模式: run-{timestamp}-{run_id}，只有3段
        assert len(parts) == 3

    def test_small_dir_max_length_truncation(self, tmp_path, monkeypatch):
        """较小的 dir_max_length 时目录名字节长度仍不超限"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", MagicMock())

        run_dir = ensure_run_dir(tmp_path, "a" * 512, retry_interval=0.01, dir_max_length=50)

        assert len(run_dir.name.encode("utf-8")) <= 50

    def test_chinese_hostname_byte_length_within_limit(self, tmp_path, monkeypatch):
        """中文 hostname 的目录名 UTF-8 字节长度不超限"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.socket.gethostname", lambda: "我的电脑" * 20)
        monkeypatch.setattr("swanlab.sdk.cmd.init.fork.current_pid", lambda: 12345)

        run_dir = ensure_run_dir(tmp_path, "run42", parallel="shared", retry_interval=0.01, dir_max_length=80)

        assert len(run_dir.name.encode("utf-8")) <= 80

    def test_chinese_run_id_byte_length_within_limit(self, tmp_path, monkeypatch):
        """包含中文的 run_id 截断后字节长度不超限"""
        monkeypatch.setattr("swanlab.sdk.cmd.init.console.warning", MagicMock())

        run_dir = ensure_run_dir(tmp_path, "测试" * 100, retry_interval=0.01, dir_max_length=80)

        assert len(run_dir.name.encode("utf-8")) <= 80
