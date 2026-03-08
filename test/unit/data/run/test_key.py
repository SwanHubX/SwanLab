"""
@author: cunyue
@file: test_key.py
@time: 2025/7/2 14:25
@description: 测试 SwanLabKey 的相关功能
"""

import pytest

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

    @pytest.mark.parametrize(
        "key, expected_section",
        [
            ("loss/metrics/label1", "loss/metrics"),  # 多个斜杠，取最后一个斜杠前的部分
            ("train/loss", "train"),                   # 单个斜杠
            ("accuracy", "default"),                   # 无斜杠，保持默认
            ("a/b/c/d", "a/b/c"),                     # 更多斜杠
        ],
    )
    def test_section_split(self, key, expected_section):
        """
        测试 section 按最后一个斜杠分割的逻辑
        """
        with UseMockRunState() as run_state:
            # 创建 DataWrapper 并解析
            data = DataWrapper(key=key, data=[Line(0.5)])
            data.parse(step=0, key=key)

            # 创建 SwanLabKey 对象
            swanlab_key = SwanLabKey(key=key, media_dir=run_state.store.media_dir, log_dir=run_state.store.log_dir)

            # 创建列信息
            column_info = swanlab_key.create_column(
                key=key,
                name=None,
                column_class="CUSTOM",
                column_config=None,
                section_type="PUBLIC",
                data=data,
                num=0
            )

            # 验证 section_name
            assert column_info.section_name == expected_section
