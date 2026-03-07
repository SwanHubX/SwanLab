"""
@author: cunyue
@file: test_apikey_helper.py
@time: 2026/3/7 14:48
@description: 测试API Key管理工具函数
"""

from swanlab.sdk.internal.apikey.helper import remove_host_suffix


class TestRemoveHostSuffix:
    """测试 remove_host_suffix 函数"""

    def test_remove_existing_suffix(self):
        """测试移除存在的后缀"""
        assert remove_host_suffix("example.com/api", "/api") == "example.com"

    def test_suffix_not_exists(self):
        """测试后缀不存在时保持不变"""
        assert remove_host_suffix("example.com", "/api") == "example.com"

    def test_empty_suffix(self):
        """测试空后缀返回原host"""
        assert remove_host_suffix("example.com", "") == "example.com"

    def test_host_with_trailing_spaces(self):
        """测试自动去除host尾部空格"""
        assert remove_host_suffix("example.com  ", "") == "example.com"

    def test_remove_suffix_with_trailing_spaces(self):
        """测试带尾部空格的host移除后缀"""
        assert remove_host_suffix("example.com/api  ", "/api") == "example.com"

    def test_partial_match_not_removed(self):
        """测试部分匹配不会被移除"""
        assert remove_host_suffix("example.com/api", "/api/v1") == "example.com/api"

    def test_remove_with_port(self):
        """测试带端口的host移除后缀"""
        assert remove_host_suffix("https://example.com:8080/api", "/api") == "https://example.com:8080"

    def test_whitespace_only_host(self):
        """测试只有空格的host"""
        assert remove_host_suffix("   ", "/api") == ""

    def test_suffix_in_middle_not_removed(self):
        """测试后缀在中间不会被移除"""
        assert remove_host_suffix("example.com/api/v1", "/api") == "example.com/api/v1"
