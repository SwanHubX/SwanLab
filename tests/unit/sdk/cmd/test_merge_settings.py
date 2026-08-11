"""
@author: cunyue
@file: test_merge_settings.py
@time: 2026/3/14
@description: 测试 swanlab.merge_settings 方法
"""

from swanlab.sdk.cmd.merge_settings import merge_settings
from swanlab.sdk.internal.settings import Settings, create_settings


class TestMergeSettings:
    def test_merge_from_dict(self):
        """传入 dict 时，后续配置快照应按键更新"""
        assert create_settings().mode != "offline"  # 确认初始状态非 offline

        merge_settings({"mode": "offline"})

        assert create_settings().mode == "offline"

    def test_merge_from_settings_object(self):
        """传入 Settings 对象时，后续配置快照应按字段更新"""
        custom = Settings(mode="local")

        merge_settings(custom)

        assert create_settings().mode == "local"

    def test_create_settings_returns_independent_snapshots(self):
        """create_settings 应返回继承全局覆盖的独立快照"""
        merge_settings({"mode": "disabled"})

        first = create_settings()
        second = create_settings()
        assert first is not second
        assert first.mode == "disabled"
        assert second.mode == "disabled"

    def test_multiple_merges_are_cumulative(self):
        """多次调用 merge_settings 应累积生效"""
        merge_settings({"mode": "offline"})
        merge_settings({"mode": "disabled"})

        assert create_settings().mode == "disabled"

    def test_api_host_derives_web_host(self):
        """仅传 api_host 时，web_host 应自动从 api_host 推导"""
        merge_settings({"api_host": "https://priv.x.cn"})

        s = create_settings()
        assert s.api_host == "https://priv.x.cn"
        assert s.web_host == "https://priv.x.cn"

    def test_api_host_with_web_host_keeps_both(self):
        """同时传 api_host 和 web_host 时，各自独立、不推导"""
        merge_settings({"api_host": "https://api.x.cn", "web_host": "https://web.x.cn"})

        s = create_settings()
        assert s.api_host == "https://api.x.cn"
        assert s.web_host == "https://web.x.cn"

    def test_repeated_api_host_merge_re_derives_web_host(self):
        """连续两次仅传 api_host 时，web_host 应每次重新推导"""
        merge_settings({"api_host": "https://a.x.cn"})
        assert create_settings().web_host == "https://a.x.cn"

        merge_settings({"api_host": "https://b.x.cn"})
        s = create_settings()
        assert s.api_host == "https://b.x.cn"
        assert s.web_host == "https://b.x.cn"
