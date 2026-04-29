"""
@author: cunyue
@file: test_echarts.py
@time: 2026/4/29
@description: ECharts TransformMedia 单元测试
"""

import json
from pathlib import Path

import pyecharts.charts
import pytest

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.sdk.internal.run.transforms.echarts import ECharts


class TestEChartsInit:
    def test_accepts_pyechnarts_chart(self):
        chart = pyecharts.charts.Bar().add_xaxis(["a", "b"]).add_yaxis("series", [1, 2])
        e = ECharts(chart)
        assert e.json_content is not None

    def test_rejects_non_chart_object(self):
        with pytest.raises(TypeError, match="Unsupported chart type"):
            ECharts("not a chart")  # type: ignore

    def test_rejects_object_without_dump_options(self):
        with pytest.raises(TypeError, match="Unsupported chart type"):
            ECharts(42)  # type: ignore

    def test_caption_stored(self):
        chart = pyecharts.charts.Line()
        e = ECharts(chart, caption="my chart")
        assert e.caption == "my chart"


class TestEChartsColumnType:
    def test_column_type_is_echarts(self):
        assert ECharts.column_type() == ColumnType.COLUMN_TYPE_ECHARTS


class TestEChartsTransform:
    def test_writes_json_file(self, tmp_path: Path):
        chart = pyecharts.charts.Bar().add_xaxis(["x"]).add_yaxis("y", [1])
        e = ECharts(chart)
        item = e.transform(step=0, path=tmp_path)
        assert (tmp_path / item.filename).exists()
        assert item.filename.endswith(".json")

    def test_json_content_valid(self, tmp_path: Path):
        chart = pyecharts.charts.Bar().add_xaxis(["x"]).add_yaxis("y", [1])
        e = ECharts(chart)
        item = e.transform(step=0, path=tmp_path)
        content = (tmp_path / item.filename).read_text()
        parsed = json.loads(content)
        # pyecharts dump_options 返回含 xAxis 等键的 dict
        assert isinstance(parsed, dict)

    def test_sha256_matches_content(self, tmp_path: Path):
        import hashlib

        chart = pyecharts.charts.Pie()
        e = ECharts(chart)
        item = e.transform(step=0, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        expected = hashlib.sha256(content).hexdigest()
        assert item.sha256 == expected

    def test_size_matches_content(self, tmp_path: Path):
        chart = pyecharts.charts.Line()
        e = ECharts(chart)
        item = e.transform(step=0, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        assert item.size == len(content)

    def test_nesting_unwrap(self):
        """ECharts 套娃：传入已有 ECharts 实例时解包内部数据"""
        chart = pyecharts.charts.Bar().add_xaxis(["x"]).add_yaxis("y", [1])
        inner = ECharts(chart)
        outer = ECharts(inner)
        assert outer.json_content == inner.json_content
