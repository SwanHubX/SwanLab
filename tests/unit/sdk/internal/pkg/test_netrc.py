"""
@author: cunyue
@file: test_netrc.py
@time: 2026/3/8 15:57
@description: 测试Netrc工具函数
"""

import netrc
import os
import sys

import pytest

from swanlab.sdk.internal.pkg.netrc import (
    get_nrc_path,
    read_netrc_by_host,
    remove_host_suffix,
    write_netrc,
)


def test_get_nrc_path(tmp_path):
    """测试获取 netrc 文件路径"""
    expected = tmp_path / ".netrc"
    assert get_nrc_path(tmp_path) == expected


def test_remove_host_suffix():
    """测试移除 host 后缀"""
    assert remove_host_suffix("api.swanlab.cn/api", "/api") == "api.swanlab.cn"
    assert remove_host_suffix("api.swanlab.cn/api/", "/api/") == "api.swanlab.cn"
    # 不匹配的后缀不作处理
    assert remove_host_suffix("api.swanlab.cn", "/api") == "api.swanlab.cn"
    # 多个后缀依次匹配（取第一个命中的）
    assert remove_host_suffix("api.swanlab.cn/v1", "/api", "/v1") == "api.swanlab.cn"
    # 空后缀处理
    assert remove_host_suffix("api.swanlab.cn", "") == "api.swanlab.cn"


class TestReadNetrc:
    def test_read_not_exists(self, tmp_path):
        """文件不存在时返回 None"""
        nrc_path = tmp_path / ".netrc"
        assert read_netrc_by_host(nrc_path, "api.swanlab.cn") is None

    def test_read_normal(self, tmp_path):
        """正常读取存在的凭证"""
        nrc_path = tmp_path / ".netrc"
        nrc_path.write_text("machine api.swanlab.cn login test_user password test_key")

        result = read_netrc_by_host(nrc_path, "api.swanlab.cn")
        assert result == ("test_user", "test_key")

        # 读取不存在的 host
        assert read_netrc_by_host(nrc_path, "other.host.com") is None

    def test_read_legacy_api_fallback(self, tmp_path):
        """测试向下兼容：读取带 /api 的旧配置，并自动清洗回写文件"""
        nrc_path = tmp_path / ".netrc"
        nrc_path.write_text("machine api.swanlab.cn/api login old_user password old_key")

        # 请求不带 /api 的纯净 host
        result = read_netrc_by_host(nrc_path, "api.swanlab.cn")

        # 1. 验证返回值正确
        assert result == ("old_user", "old_key")

        # 2. 验证文件已经被自愈修复（移除了 /api，且强覆盖）
        parsed = netrc.netrc(nrc_path)
        assert "api.swanlab.cn/api" not in parsed.hosts
        assert "api.swanlab.cn" in parsed.hosts
        assert parsed.authenticators("api.swanlab.cn") == ("old_user", "", "old_key")

    def test_read_exceptions(self, tmp_path):
        """测试读取时的异常处理（目录冲突和格式损坏）"""
        nrc_path = tmp_path / ".netrc"

        # 1. 格式损坏
        nrc_path.write_text("invalid content without machine keyword")
        with pytest.raises(IOError, match="Failed to access or parse netrc file"):
            read_netrc_by_host(nrc_path, "api.swanlab.cn")

        # 2. 目录冲突
        nrc_path.unlink()
        nrc_path.mkdir()
        with pytest.raises(IOError, match="Failed to access or parse netrc file"):
            read_netrc_by_host(nrc_path, "api.swanlab.cn")


class TestWriteNetrc:
    def test_write_normal_and_permissions(self, tmp_path):
        """测试正常写入和权限控制"""
        nrc_path = tmp_path / ".netrc"
        write_netrc(nrc_path, "api.swanlab.cn", "test_user", "test_key")

        assert nrc_path.exists()

        # 权限校验：仅在非 Windows 系统上执行严格比对 (0o600)
        if sys.platform != "win32":
            assert (nrc_path.stat().st_mode & 0o777) == 0o600
        else:
            assert os.access(nrc_path, os.R_OK | os.W_OK)

        # 验证写入内容
        parsed = netrc.netrc(nrc_path)
        assert parsed.authenticators("api.swanlab.cn") == ("test_user", "", "test_key")

    def test_write_global_sso(self, tmp_path):
        """测试全局单点登录机制：后写入的必须清空之前不同环境的所有凭证"""
        nrc_path = tmp_path / ".netrc"

        # 第一次写入（公有云）
        write_netrc(nrc_path, "api.swanlab.cn", "user1", "key1")

        # 第二次写入（私有化节点），验证是否顶掉了公有云的凭证
        write_netrc(nrc_path, "private.company.com", "user2", "key2")

        parsed = netrc.netrc(nrc_path)

        # 验证新节点已生效
        assert parsed.authenticators("private.company.com") == ("user2", "", "key2")
        # 验证旧节点已被强行清除（只保留了 1 个 host）
        assert "api.swanlab.cn" not in parsed.hosts
        assert len(parsed.hosts) == 1

    def test_write_recovery(self, tmp_path):
        """测试当 netrc 文件内容损坏时，写入操作能够重置并成功保存"""
        nrc_path = tmp_path / ".netrc"
        nrc_path.write_text("broken file contents!!")

        # 即使文件损坏，写入逻辑也应当能捕获异常并覆盖
        write_netrc(nrc_path, "api.swanlab.cn", "new_user", "new_key")

        parsed = netrc.netrc(nrc_path)
        assert parsed.authenticators("api.swanlab.cn") == ("new_user", "", "new_key")

    def test_write_directory_conflict(self, tmp_path):
        """测试目标路径被目录占用时的防御策略"""
        nrc_path = tmp_path / ".netrc"
        nrc_path.mkdir()

        with pytest.raises(IOError, match="is a directory, but a file is expected"):
            write_netrc(nrc_path, "api.swanlab.cn", "user", "key")
