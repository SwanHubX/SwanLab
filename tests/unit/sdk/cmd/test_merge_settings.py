"""
@author: cunyue
@file: test_merge_settings.py
@time: 2026/3/14
@description: 测试 swanlab.merge_settings 方法
"""

from swanlab.sdk.cmd.merge_settings import merge_settings
from swanlab.sdk.internal.settings import Settings, settings


class TestMergeSettings:
    def test_merge_from_dict(self):
        """传入 dict 时，全局 settings 应按键更新"""
        assert settings.mode != "offline"  # 确认初始状态非 offline

        merge_settings({"mode": "offline"})

        assert settings.mode == "offline"

    def test_merge_from_settings_object(self):
        """传入 Settings 对象时，全局 settings 应按字段更新"""
        custom = Settings(mode="local")

        merge_settings(custom)

        assert settings.mode == "local"

    def test_merge_updates_global_singleton(self):
        """merge_settings 应直接修改全局单例，而不是创建副本"""
        original_id = id(settings)

        merge_settings({"mode": "disabled"})

        # 对象身份不变（仍是同一个 settings 实例）
        assert id(settings) == original_id
        assert settings.mode == "disabled"

    def test_multiple_merges_are_cumulative(self):
        """多次调用 merge_settings 应累积生效"""
        merge_settings({"mode": "offline"})
        merge_settings({"mode": "disabled"})

        assert settings.mode == "disabled"
