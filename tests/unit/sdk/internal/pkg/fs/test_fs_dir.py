"""
@author: cunyue
@file: test_fs_dir.py
@time: 2026/3/11 14:30
@description: SwanLab SDK 文件系统辅助函数测试
"""

import importlib
import sys
from pathlib import Path

import pytest

from swanlab.sdk.internal.pkg.fs import dir

sys.modules[".."] = type(sys)("mock_parent")
sys.modules[".."].console = type(sys)("mock_console")


def test_safe_mkdir_normal(tmp_path: Path):
    """测试常规单级目录创建"""
    target = tmp_path / "test_dir"
    result = dir.safe_mkdir(target)

    assert result == target
    assert target.exists()
    assert target.is_dir()


def test_safe_mkdir_nested(tmp_path: Path):
    """测试多级目录创建"""
    target = tmp_path / "a" / "b" / "c"
    dir.safe_mkdir(target)

    assert target.exists()
    assert target.is_dir()


def test_safe_mkdir_already_exists(tmp_path: Path):
    """测试目录已存在时的情况"""
    target = tmp_path / "exist_dir"
    target.mkdir()

    # 应该直接通过，不会报错
    result = dir.safe_mkdir(target)
    assert result == target


def test_safe_mkdir_nas_timeout(monkeypatch, tmp_path: Path):
    """测试极端 NAS 延迟/权限问题下，探针是否能正确抛出 TimeoutError"""
    target = tmp_path / "timeout_dir"

    time_calls = [0, 10, 20, 30]
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.time.time", lambda: time_calls.pop(0))

    def mock_tempfile(*args, **kwargs):
        raise OSError("Simulated NAS Permission Denied")

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.tempfile.TemporaryFile", mock_tempfile)

    with pytest.raises(TimeoutError, match="is not writable within"):
        dir.safe_mkdir(target, timeout=5.0)


def test_timeout_env_invalid_string(monkeypatch):
    """测试环境变量为非法字符串时，是否能安全回退到默认值 5.0"""
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "invalid_abc123")

    mock_warnings = []
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.console.warning", lambda msg: mock_warnings.append(msg))

    importlib.reload(dir)

    assert dir.TIMEOUT == 5.0
    assert len(mock_warnings) == 1
    assert "Invalid SWANLAB_FS_TIMEOUT" in mock_warnings[0]

    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)


def test_timeout_env_negative_number(monkeypatch):
    """测试环境变量为负数或0时，是否能安全回退到默认值 5.0"""
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "-10.5")

    mock_warnings = []
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.console.warning", lambda msg: mock_warnings.append(msg))

    importlib.reload(dir)

    assert dir.TIMEOUT == 5.0
    assert len(mock_warnings) == 1
    assert "must be > 0" in mock_warnings[0]

    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)
