"""
@author: cunyue
@file: test_apikey.py
@time: 2026/3/7 15:24
@description: 测试API Key管理
"""

import netrc
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal import apikey


@pytest.fixture
def mock_env(mocker, tmp_path):
    """
    统一准备测试环境的 fixture，拦截外部依赖，防止污染本地真实配置
    """
    # 1. Mock Settings
    mock_settings = MagicMock()
    mock_settings.api_url = "https://api.swanlab.cn/api"
    mock_settings.api_key = None  # 默认不设置，方便测试文件读取逻辑
    mock_settings.save_dir = tmp_path / "swanlab"

    # 替换为真实的模块路径
    mocker.patch("swanlab.sdk.internal.apikey.get_current_settings", return_value=mock_settings)

    # 2. Mock nrc 文件路径，指向 pytest 提供的安全的临时文件
    nrc_file = tmp_path / ".netrc"
    mocker.patch("swanlab.sdk.internal.apikey.get_nrc_path", return_value=nrc_file)

    # 3. Mock 字符串处理辅助函数
    mocker.patch("swanlab.sdk.internal.apikey.remove_host_suffix", return_value="api.swanlab.cn")

    return mock_settings, nrc_file


def test_get_from_settings(mock_env):
    """测试优先级最高的情况：Settings 中已存在 API Key"""
    mock_settings, _ = mock_env
    mock_settings.api_key = "test-key-from-settings"
    assert apikey.exists() is True
    assert apikey.get() == "test-key-from-settings"


def test_save_and_get_from_netrc(mock_env):
    """测试完整流程：保存 API Key 到文件，并从中读取"""
    mock_settings, nrc_file = mock_env

    # 模拟用户调用保存
    apikey.save(username="test_user", api_key="test-key-from-file")

    # 验证文件是否正确生成并包含了正确的格式
    assert nrc_file.exists()
    parsed_nrc = netrc.netrc(nrc_file)
    assert parsed_nrc.authenticators("api.swanlab.cn") == ("test_user", "", "test-key-from-file")

    # 验证内存中的 settings 被更新
    mock_settings.merge_settings.assert_called_once_with({"api_key": "test-key-from-file"})

    # 验证 get() 能够成功从文件中读出来
    assert apikey.exists() is True
    assert apikey.get() == "test-key-from-file"


def test_netrc_backward_compatibility(mock_env):
    """测试向下兼容：当 netrc 中的 host 带有 /api 后缀时，能够正确读取并修复文件"""
    _, nrc_file = mock_env

    # 手动构造一个旧版本格式的 netrc 文件
    nrc_file.write_text("machine api.swanlab.cn/api login old_user password old_key")

    # 调用 get() 应当能读取成功
    assert apikey.get() == "old_key"

    # 断言文件已经被自动修复：后缀被移除
    parsed_nrc = netrc.netrc(nrc_file)
    auth = parsed_nrc.authenticators("api.swanlab.cn")
    assert auth is not None
    # 旧版本的 netrc 文件中，auth[1] 可能是 None 而非空字符串，并且我们并不关心 auth[1] 的值，所以这里不比较即可
    # assert parsed_nrc.authenticators("api.swanlab.cn") == ("old_user", "", "old_key")
    assert (auth[0], auth[2]) == ("old_user", "old_key")


def test_get_not_found(mock_env):
    """测试所有方式都找不到 API Key 的情况"""
    _, _ = mock_env
    # 默认 mock_env 既没有 settings 也没有文件
    assert apikey.exists() is False

    with pytest.raises(FileNotFoundError, match="The API key file or target host does not exist"):
        apikey.get()
