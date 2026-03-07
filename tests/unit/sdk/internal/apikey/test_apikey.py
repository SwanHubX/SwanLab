"""
@author: cunyue
@file: test_apikey.py
@time: 2026/3/7 15:24
@description: 测试API Key管理
"""

import netrc
import os
import sys
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal import apikey


@pytest.fixture
def mock_env(mocker, tmp_path, monkeypatch):
    """
    统一准备测试环境，拦截外部依赖
    """
    monkeypatch.setenv("HOME", str(tmp_path))
    mock_settings = MagicMock(api_url="https://api.swanlab.cn/api", api_key=None)
    mocker.patch("swanlab.sdk.internal.apikey.get_current_settings", return_value=mock_settings)

    nrc_file = tmp_path / ".netrc"
    mocker.patch("swanlab.sdk.internal.apikey.get_nrc_path", return_value=nrc_file)
    mocker.patch("swanlab.sdk.internal.apikey.remove_host_suffix", return_value="api.swanlab.cn")

    return mock_settings, nrc_file


def test_save_and_get_flow(mock_env):
    """测试完整流程：保存、权限校验、解构读取"""
    mock_settings, nrc_file = mock_env
    apikey.save(username="test_user", api_key="test-key-from-file")

    # 验证文件权限与内容
    assert nrc_file.exists()
    # 权限校验：仅在非 Windows 系统上执行数值比对
    if sys.platform != "win32":
        assert (nrc_file.stat().st_mode & 0o777) == 0o600
    else:
        # Windows 上仅验证文件是否仍可读
        assert os.access(nrc_file, os.R_OK)

    parsed_nrc = netrc.netrc(nrc_file)
    auth = parsed_nrc.authenticators("api.swanlab.cn")
    assert auth is not None
    # 旧版本的 netrc 文件中，auth[1] 可能是 None 而非空字符串，并且我们并不关心 auth[1] 的值，所以这里不比较即可
    # assert parsed_nrc.authenticators("api.swanlab.cn") == ("test_user", "", "test-key-from-file")
    assert (auth[0], auth[2]) == ("test_user", "test-key-from-file")

    # 验证内存更新与获取
    mock_settings.merge_settings.assert_called_once_with({"api_key": "test-key-from-file"})
    assert apikey.get() == "test-key-from-file"


def test_directory_conflict(mock_env):
    """测试路径被目录占用时的防御"""
    _, nrc_file = mock_env
    nrc_file.mkdir()  # 创建同名目录
    with pytest.raises(IOError, match="is a directory"):
        apikey.save("user", "key")


def test_corrupted_file_recovery(mock_env):
    """测试 netrc 文件损坏时能自动重置并保存"""
    _, nrc_file = mock_env
    nrc_file.write_text("invalid netrc content...")
    apikey.save("new_user", "new_key")
    assert apikey.get() == "new_key"


def test_netrc_backward_compatibility(mock_env):
    """测试旧版 host (/api) 的读取与自动修复"""
    _, nrc_file = mock_env
    nrc_file.write_text("machine api.swanlab.cn/api login old_user password old_key")

    assert apikey.get() == "old_key"
    # 验证修复：后缀应已被移除
    auth = netrc.netrc(nrc_file).authenticators("api.swanlab.cn")
    assert auth is not None
    assert (auth[0], auth[2]) == ("old_user", "old_key")


def test_get_not_found(mock_env):
    """测试找不到 API Key 时的异常"""
    assert apikey.exists() is False
    with pytest.raises(FileNotFoundError, match="does not exist"):
        apikey.get()
