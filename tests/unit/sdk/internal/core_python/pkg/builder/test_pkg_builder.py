"""
@author: cunyue
@file: test_pkg_builder.py
@time: 2026/5/8
@description: 测试 swanlab.sdk.internal.pkg.builder 模块
"""

import pytest

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnType, SectionType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.sdk.internal.core_python.context import CoreConfig, CoreContext
from swanlab.sdk.internal.core_python.pkg import builder


def make_ctx(tmp_path, section_rule: int) -> CoreContext:
    return CoreContext(
        config=CoreConfig(
            run_id="test-run-id",
            run_dir=tmp_path / "run",
            section_rule=section_rule,
            record_batch=10000,
            record_interval=5.0,
            save_split=100 * 1024 * 1024,
            save_size=50 * 1024 * 1024 * 1024,
            save_part=32 * 1024 * 1024,
            save_batch=100,
        )
    )


class TestBuildAutoColumn:
    @pytest.mark.parametrize(
        "index, expected_section",
        [
            (0, "train"),
            (-1, "train/loss"),
            (5, "train/loss"),
            (6, "train"),
        ],
    )
    def test_scalar_section_rule(self, tmp_path, index, expected_section):
        ctx = make_ctx(tmp_path, index)
        record = ScalarRecord(key="train/loss/epoch", type=ColumnType.COLUMN_TYPE_SCALAR)

        column = builder.build_auto_column(ctx, record)

        assert column.column_class == ColumnClass.COLUMN_CLASS_CUSTOM
        assert column.column_key == "train/loss/epoch"
        assert column.column_type == ColumnType.COLUMN_TYPE_SCALAR
        assert column.section_name == expected_section
        assert column.section_type == SectionType.SECTION_TYPE_PUBLIC
        assert column.metric_name == ""
        assert column.chart_name == ""

    def test_media_section_rule(self, tmp_path):
        ctx = make_ctx(tmp_path, -1)
        record = MediaRecord(key="train/image/epoch", type=ColumnType.COLUMN_TYPE_IMAGE)

        column = builder.build_auto_column(ctx, record)

        assert column.column_key == "train/image/epoch"
        assert column.column_type == ColumnType.COLUMN_TYPE_IMAGE
        assert column.section_name == "train/image"

    def test_key_without_slash(self, tmp_path):
        ctx = make_ctx(tmp_path, -1)
        record = ScalarRecord(key="loss", type=ColumnType.COLUMN_TYPE_SCALAR)

        column = builder.build_auto_column(ctx, record)

        assert column.section_name == ""


@pytest.mark.parametrize("system", [False, True])
def test_build_resume_histogram_column(system):
    column = builder.build_resume_column("weights", system=system, column_type=ColumnType.COLUMN_TYPE_HISTOGRAM)

    assert column.column_key == "weights"
    assert column.column_type == ColumnType.COLUMN_TYPE_HISTOGRAM
    expected_class = ColumnClass.COLUMN_CLASS_SYSTEM if system else ColumnClass.COLUMN_CLASS_CUSTOM
    assert column.column_class == expected_class


def test_build_resume_scalar_and_media_keep_existing_types():
    assert builder.build_resume_column("loss").column_type == ColumnType.COLUMN_TYPE_SCALAR
    assert builder.build_resume_column("image", media=True).column_type == ColumnType.COLUMN_TYPE_UNSPECIFIED
