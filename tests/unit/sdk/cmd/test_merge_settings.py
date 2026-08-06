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
