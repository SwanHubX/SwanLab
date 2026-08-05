"""
@author: cunyue
@file: test_settings_experiment.py
@time: 2026/3/6 12:26
@description: 测试 SwanLab 实验配置
"""

import pytest
from pydantic import ValidationError

from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings.experiment import map_resume_value

QUOTE_NORMALIZATION_CASES = [
    ("SWANLAB_PROJECT_NAME", "project.name", "quoted-project"),
    ("SWANLAB_PROJ_NAME", "project.name", "legacy-project"),
    ("SWANLAB_PROJECT_WORKSPACE", "project.workspace", "quoted_workspace"),
    ("SWANLAB_WORKSPACE", "project.workspace", "legacy_workspace"),
    ("SWANLAB_EXPERIMENT_COLOR", "experiment.color", "#123ABC"),
    ("SWANLAB_EXP_COLOR", "experiment.color", "#ABC123"),
    ("SWANLAB_RUN_ID", "run.id", "quoted-run"),
    ("SWANLAB_RUN_DIR", "run.dir", "quoted-dir"),
]


def _get_nested_attr(value, path):
    for part in path.split("."):
        value = getattr(value, part)
    return value


@pytest.mark.parametrize("env_name,path,expected", QUOTE_NORMALIZATION_CASES)
def test_marked_env_strips_double_outer_quote_pair(monkeypatch, env_name, path, expected):
    for name, _, _ in QUOTE_NORMALIZATION_CASES:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(env_name, f'"{expected}"')

    assert _get_nested_attr(Settings(), path) == expected


@pytest.mark.parametrize(
    "env_name,path,expected",
    QUOTE_NORMALIZATION_CASES[:-2],
)
def test_marked_env_strips_single_outer_quote_pair(monkeypatch, env_name, path, expected):
    for name, _, _ in QUOTE_NORMALIZATION_CASES:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(env_name, f"'{expected}'")

    assert _get_nested_attr(Settings(), path) == expected


@pytest.mark.parametrize("env_name,path,expected", QUOTE_NORMALIZATION_CASES[-2:])
def test_run_values_preserve_allowed_single_quotes(monkeypatch, env_name, path, expected):
    for name, _, _ in QUOTE_NORMALIZATION_CASES:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv(env_name, f"'{expected}'")

    assert _get_nested_attr(Settings(), path) == f"'{expected}'"


def test_strict_env_strips_only_one_outer_quote_pair(monkeypatch):
    monkeypatch.setenv("SWANLAB_RUN_ID", '""abc""')

    with pytest.raises(ValidationError):
        Settings()


def test_non_strict_env_preserves_outer_quotes(monkeypatch):
    value = "'quoted description'"
    monkeypatch.setenv("SWANLAB_EXPERIMENT_DESCRIPTION", value)

    assert Settings().experiment.description == value


def test_outer_quote_normalization_only_applies_to_env():
    value = "'quoted-run'"

    assert Settings(run={"id": value}).run.id == value


def test_experiment_tags_parse(monkeypatch):
    """
    测试实验配置中的 tags 字段解析，同时支持逗号分隔和JSON字符串
    JSON字符串必须是数组格式，如 `["tag1", "tag2"]`
    """
    # 1. 逗号分隔
    monkeypatch.setenv("SWANLAB_EXPERIMENT_TAGS", "tag1, tag2, tag3")
    s = Settings()
    assert s.experiment.tags == ["tag1", "tag2", "tag3"]
    # 2. JSON数组
    monkeypatch.setenv("SWANLAB_EXPERIMENT_TAGS", '["tag1", "tag2", "tag3"]')
    s = Settings()
    assert s.experiment.tags == ["tag1", "tag2", "tag3"]
    # 3. 空值
    monkeypatch.setenv("SWANLAB_EXPERIMENT_TAGS", "")
    s = Settings()
    assert s.experiment.tags == []
    # 4. JSON 字典
    monkeypatch.setenv("SWANLAB_EXPERIMENT_TAGS", '{"key": "value"}')
    s = Settings()
    assert s.experiment.tags == ['{"key": "value"}']
    # 5. model_validator 解析
    s = Settings.model_validate({"experiment": {"tags": ["tag1", "tag2"]}})
    assert s.experiment.tags == ["tag1", "tag2"]
    # 6. model_validator 异常
    with pytest.raises(ValueError, match="tags must be a list, dict, or string"):
        Settings.model_validate({"experiment": {"tags": 123}})


def test_run_resume_parse(monkeypatch):
    """
    测试 RunSettings.resume 字段的解析
    """
    monkeypatch.setenv("SWANLAB_RUN_RESUME", "allow")
    s = Settings()
    assert s.run.resume == "allow"
    monkeypatch.setenv("SWANLAB_RUN_RESUME", "never")
    s = Settings()
    assert s.run.resume == "never"
    monkeypatch.delenv("SWANLAB_RUN_RESUME")
    s = Settings()
    assert s.run.resume == "never"
    s = Settings.model_validate({"run": {"resume": "allow"}})
    assert s.run.resume == "allow"
    s = Settings.model_validate({"run": {"resume": "1"}})
    assert s.run.resume == "allow"
    with pytest.raises(ValueError, match="Invalid resume value"):
        Settings.model_validate({"run": {"resume": "invalid"}})
    with pytest.raises(ValueError, match="Invalid resume value"):
        Settings.model_validate({"run": {"resume": 1}})


def test_run_parallel_forces_resume_allow():
    s = Settings.model_validate({"run": {"parallel": "shared"}})
    assert s.run.resume == "allow"

    s = Settings.model_validate({"run": {"parallel": "shared", "resume": "never"}})
    assert s.run.resume == "allow"


class TestMapResumeValue:
    @pytest.mark.parametrize(
        "input_value,expected",
        [
            # Boolean values
            (True, "allow"),
            (False, "never"),
            # String literals
            ("must", "must"),
            ("allow", "allow"),
            ("never", "never"),
            # Boolean-like strings
            ("true", "allow"),
            ("false", "never"),
            ("True", "allow"),
            ("False", "never"),
            ("yes", "allow"),
            ("no", "never"),
            ("1", "allow"),
            ("0", "never"),
        ],
    )
    def test_map_resume_value(self, input_value, expected):
        """测试 map_resume_value 函数的各种输入映射"""
        assert map_resume_value(input_value) == expected

    def test_map_resume_value_invalid(self):
        """测试 map_resume_value 函数的异常情况"""
        with pytest.raises(ValueError, match="Invalid resume value"):
            map_resume_value("invalid")
        with pytest.raises(ValueError, match="Invalid resume value"):
            map_resume_value("maybe")


class TestRunDir:
    def test_dir_default_none(self):
        s = Settings()
        assert s.run.dir is None

    def test_dir_exceeds_dir_max_length(self):
        with pytest.raises(Exception, match="run.dir length"):
            Settings.model_validate({"run": {"dir": "a" * 51, "dir_max_length": 50}})

    def test_dir_within_dir_max_length(self):
        s = Settings.model_validate({"run": {"dir": "a" * 50, "dir_max_length": 50}})
        assert s.run.dir is not None
        assert len(s.run.dir) == 50

    def test_dir_env_injection(self, monkeypatch):
        monkeypatch.setenv("SWANLAB_RUN_DIR", "from-env")
        s = Settings()
        assert s.run.dir == "from-env"
