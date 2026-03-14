"""
@author: cunyue
@file: test_key.py
@time: 2025/7/2 14:25
@description: 测试 SwanLabKey 的相关功能
"""

import pytest

from swanlab.data.modules import DataWrapper, Line
from swanlab.data.run.key import SwanLabKey
from swanlab.data.modules import Line, DataWrapper
from swanlab.toolkit import ChartType
from tutils.setup import UseMockRunState


class TestMockKey:
    """
    测试模拟一个 SwanLabKey 对象
    """

    @pytest.mark.parametrize(
        "column_class, expected_section_type",
        [
            ["CUSTOM", "PUBLIC"],
            ["SYSTEM", "SYSTEM"],
        ],
    )
    def test_column_class_ok(self, column_class, expected_section_type):
        """
        测试不同的 column_class 是否能正确映射到预期的 SectionType
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type="FLOAT",
                column_class=column_class,
                error=None,
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=None,
            )
            assert column_info.section_type == expected_section_type

    @pytest.mark.parametrize(
        "column_type, expected_section_type",
        [
            ["FLOAT", "PUBLIC"],
            ["IMAGE", "PUBLIC"],
            ["AUDIO", "PUBLIC"],
            ["TEXT", "PUBLIC"],
            ["OBJECT3D", "PUBLIC"],
            ["MOLECULE", "PUBLIC"],
            ["ECHARTS", "CUSTOM"],
        ],
    )
    def test_column_type_section_type(self, column_type, expected_section_type):
        """
        测试不同的 column_type 是否能正确映射到预期的 SectionType
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type=column_type,
                column_class="CUSTOM",
                error=None,
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=None,
            )
            assert column_info.section_type == expected_section_type

    @pytest.mark.parametrize(
        "column_type, expected_chart_type",
        [
            ["FLOAT", ChartType.LINE],
            ["IMAGE", ChartType.IMAGE],
            ["AUDIO", ChartType.AUDIO],
            ["TEXT", ChartType.TEXT],
            ["OBJECT3D", ChartType.OBJECT3D],
            ["MOLECULE", ChartType.MOLECULE],
            ["ECHARTS", ChartType.ECHARTS],
        ],
    )
    def test_column_type_ok(self, column_type, expected_chart_type):
        """
        测试不同的 column_type 是否能正确映射到预期的图表类型
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type=column_type,
                column_class="CUSTOM",
                error=None,
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=None,
            )
            assert column_info.chart_type == expected_chart_type

    def test_column_unknown_type(self):
        """
        未知的 column_type 应该报错
        """
        with UseMockRunState() as run_state:
            with pytest.raises(RuntimeError) as exc_info:
                SwanLabKey.mock_from_remote(
                    key="test",
                    column_type="UNKNOWN",
                    column_class="CUSTOM",
                    error=None,
                    media_dir=run_state.store.media_dir,
                    log_dir=run_state.store.log_dir,
                    kid=0,
                    step=None,
                )
            assert (
                str(exc_info.value)
                == "Unknown chart type: UNKNOWN, maybe you need to update swanlab: pip install -U swanlab"
            )

    def test_column_with_error(self):
        """
        测试传递列的错误信息时是否能正确处理
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type="FLOAT",
                column_class="CUSTOM",
                error={"excepted": "float", "data_class": "string"},
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=None,
            )
            assert key_obj.column_info.error.expected == "float"
            assert key_obj.column_info.error.got == "string"

    def test_column_with_error_unknown(self):
        """
        测试传递错误信息时，如果 excepted 或 data_class 为 None，应该报错
        因为这不符合目前的错误格式要求
        """
        with UseMockRunState() as run_state:
            with pytest.raises(RuntimeError) as exc_info:
                SwanLabKey.mock_from_remote(
                    key="test",
                    column_type="FLOAT",
                    column_class="CUSTOM",
                    error={"excepted": None, "data_class": "string"},
                    media_dir=run_state.store.media_dir,
                    log_dir=run_state.store.log_dir,
                    kid=0,
                    step=None,
                )
            assert (
                str(exc_info.value)
                == "Invalid error format: {'excepted': None, 'data_class': 'string'}, expected and got must be provided. "
                "Maybe you need to update swanlab: pip install -U swanlab"
            )

    def test_step_none(self):
        """
        测试 step 为 None 的情况
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type="FLOAT",
                column_class="CUSTOM",
                error=None,
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=None,
            )
            assert len(key_obj.steps) == 0, "Steps should be empty when step is None"

    def test_step_not_none(self):
        """
        测试 step 不为 None 的情况
        """
        with UseMockRunState() as run_state:
            key_obj, column_info = SwanLabKey.mock_from_remote(
                key="test",
                column_type="FLOAT",
                column_class="CUSTOM",
                error=None,
                media_dir=run_state.store.media_dir,
                log_dir=run_state.store.log_dir,
                kid=0,
                step=100,
            )
            assert len(key_obj.steps) == 101, "Steps should contain one entry when step is provided"
            assert 1 in key_obj.steps, "Step 1 should be present in the steps"
            assert 101 not in key_obj.steps, "Step 10 should be present in the steps"


class TestKeySummary:
    @staticmethod
    def _new_key(run_state, key: str = "test") -> SwanLabKey:
        return SwanLabKey(key, run_state.store.media_dir, run_state.store.log_dir)

    @staticmethod
    def _add_line(key_obj: SwanLabKey, value, step: int):
        data = DataWrapper(key_obj.key, [Line(value)])
        data.parse(step=step, key=key_obj.key)
        if not key_obj.chart_created:
            key_obj.create_column(
                key=key_obj.key,
                name=None,
                column_class="CUSTOM",
                column_config=None,
                section_type="PUBLIC",
                data=data,
                num=0,
            )
        return key_obj.add(data)

    def test_overwrite_non_extreme_promote_to_max_without_rebuild(self, monkeypatch):
        with UseMockRunState() as run_state:
            key_obj = self._new_key(run_state)
            self._add_line(key_obj, 1, 0)
            self._add_line(key_obj, 10, 1)
            self._add_line(key_obj, 5, 2)

            monkeypatch.setattr(key_obj, "_rebuild_summary", lambda: pytest.fail("unexpected summary rebuild"))
            metric_info = self._add_line(key_obj, 20, 2)

            assert metric_info.metric_overwrite is True
            assert metric_info.metric_summary == {
                "max": 20.0,
                "max_step": 2,
                "min": 1.0,
                "min_step": 0,
                "num": 3,
            }

    def test_overwrite_current_max_demote_triggers_rebuild(self, monkeypatch):
        with UseMockRunState() as run_state:
            key_obj = self._new_key(run_state)
            self._add_line(key_obj, 1, 0)
            self._add_line(key_obj, 10, 1)
            self._add_line(key_obj, 5, 2)

            rebuild_calls = []
            original_rebuild = key_obj._rebuild_summary

            def rebuild():
                rebuild_calls.append(True)
                original_rebuild()

            monkeypatch.setattr(key_obj, "_rebuild_summary", rebuild)
            metric_info = self._add_line(key_obj, 7, 1)

            assert len(rebuild_calls) == 1
            assert metric_info.metric_summary == {
                "max": 7.0,
                "max_step": 1,
                "min": 1.0,
                "min_step": 0,
                "num": 3,
            }

    def test_overwrite_duplicate_extreme_non_owner_skips_rebuild(self, monkeypatch):
        with UseMockRunState() as run_state:
            key_obj = self._new_key(run_state)
            self._add_line(key_obj, 1, 0)
            self._add_line(key_obj, 10, 1)
            self._add_line(key_obj, 10, 2)

            monkeypatch.setattr(key_obj, "_rebuild_summary", lambda: pytest.fail("duplicate extremum should not rebuild"))
            metric_info = self._add_line(key_obj, 9, 2)

            assert metric_info.metric_summary == {
                "max": 10.0,
                "max_step": 1,
                "min": 1.0,
                "min_step": 0,
                "num": 3,
            }

    def test_overwrite_current_max_to_nan_triggers_rebuild(self, monkeypatch):
        with UseMockRunState() as run_state:
            key_obj = self._new_key(run_state)
            self._add_line(key_obj, 1, 0)
            self._add_line(key_obj, 10, 1)
            self._add_line(key_obj, 5, 2)

            rebuild_calls = []
            original_rebuild = key_obj._rebuild_summary

            def rebuild():
                rebuild_calls.append(True)
                original_rebuild()

            monkeypatch.setattr(key_obj, "_rebuild_summary", rebuild)
            metric_info = self._add_line(key_obj, float("nan"), 1)

            assert len(rebuild_calls) == 1
            assert metric_info.metric_summary == {
                "max": 5.0,
                "max_step": 2,
                "min": 1.0,
                "min_step": 0,
                "num": 3,
            }
