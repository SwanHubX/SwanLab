"""
@author: cunyue
@file: test_fs_dir.py
@time: 2026/3/11 14:30
@description: SwanLab SDK 文件系统辅助函数测试
"""

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.pkg.fs import dir

sys.modules[".."] = MagicMock()
sys.modules[".."].console = MagicMock()


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


@patch("swanlab.sdk.internal.pkg.fs.dir.tempfile.TemporaryFile")
@patch("swanlab.sdk.internal.pkg.fs.dir.time.time")
def test_safe_mkdir_nas_timeout(mock_time, mock_tempfile, tmp_path: Path):
    """
    测试极端 NAS 延迟/权限问题下，探针是否能正确抛出 TimeoutError
    """
    target = tmp_path / "timeout_dir"

    # 模拟时间流逝：第一次获取当前时间为 0，之后每次获取都增加 10 秒（超过 timeout 5.0）
    mock_time.side_effect = [0, 10, 20, 30]

    # 强制让探针（TemporaryFile）每次都因为权限问题抛错
    mock_tempfile.side_effect = OSError("Simulated NAS Permission Denied")

    with pytest.raises(TimeoutError, match="is not writable within"):
        dir.safe_mkdir(target, timeout=5.0)


def test_timeout_env_invalid_string(monkeypatch):
    """
    测试环境变量为非法字符串时，是否能安全回退到默认值 5.0
    并且不会发生 ValueError 导致模块导入崩溃
    """
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "invalid_abc123")

    # 模拟 .. console 以便能够断言警告输出
    with patch("swanlab.sdk.internal.pkg.fs.dir.console.warning") as mock_warn:
        importlib.reload(dir)

        # 1. 确认回退成功
        assert dir.TIMEOUT == 5.0
        # 2. 确认触发了警告信息
        mock_warn.assert_called_once()
        assert "Invalid SWANLAB_FS_TIMEOUT" in mock_warn.call_args[0][0]

    # 清理环境
    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)


def test_timeout_env_negative_number(monkeypatch):
    """
    测试环境变量为负数或0时，是否能安全回退到默认值 5.0
    """
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "-10.5")

    with patch("swanlab.sdk.internal.pkg.fs.dir.console.warning") as mock_warn:
        importlib.reload(dir)

        # 1. 确认回退成功（不允许负数超时）
        assert dir.TIMEOUT == 5.0
        # 2. 确认触发了警告信息
        mock_warn.assert_called_once()
        assert "must be > 0" in mock_warn.call_args[0][0]

    # 清理环境
    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)
