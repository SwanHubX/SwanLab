"""
@author: cunyue
@file: test_toolkit_key.py
@time: 2025/10/10 22:26
@description: 测试工具中的
"""

from swanlab.toolkit.models import metric as M
from swanlab.toolkit.models.data import ChartType


def test_column_config():
    c = M.ColumnConfig(y_range=(0, 100), chart_name="CPU Utilization (%)", metric_name=None)
    assert c.y_range == (0, 100)
    assert c.chart_name == "CPU Utilization (%)"
    assert c.metric_name is None

    clone = c.clone(y_range=None, metric_name="12345")
    assert clone.y_range == (0, 100)
    assert clone.chart_name == c.chart_name
    assert clone.metric_name == "12345"
    clone = c.clone(y_range=(0, 50), chart_name="CPU Utilization (%)", metric_name="12345")
    assert clone.y_range == (0, 50)
    assert clone.chart_name == c.chart_name
    assert clone.metric_name == "12345"
    c = M.ColumnConfig(y_range=(0, 100), chart_name="CPU Utilization (%)", metric_name="12345", chart_index="1")
    assert c.y_range == (0, 100)
    assert c.chart_name == "CPU Utilization (%)"
    assert c.metric_name == "12345"
    assert c.chart_index == "1"
    c = c.clone(metric_color=("red", "blue"))
    assert c.y_range == (0, 100)
    assert c.chart_name == "CPU Utilization (%)"
    assert c.metric_name == "12345"
    assert c.chart_index == "1"
    assert c.metric_color == ("red", "blue")


def test_column_info():
    c = M.ColumnInfo(
        key="a/1",
        kid="b",
        name="c",
        cls="SYSTEM",
        section_name="e",
        section_sort=1,
        section_type="PUBLIC",
        chart_type=ChartType.TEXT,
        chart_reference="STEP",
        error=None,
        config=None,
    )
    assert c.got is None
    assert c.key == "a/1"
    assert c.kid == "b"
    assert c.name == "c"
    assert c.cls == "SYSTEM"
    assert c.section_name == "e"
    assert c.section_sort == 1
    assert c.chart_type == ChartType.TEXT
    assert c.chart_reference == "STEP"
    assert c.section_type == "PUBLIC"
    assert c.error is None
    assert c.config is None
    assert c.key_encode == "a%2F1"


def test_metric_info():
    c = M.ColumnInfo(
        key="a/1",
        kid="b",
        name="c",
        cls="SYSTEM",
        section_name="e",
        section_sort=1,
        chart_type=ChartType.TEXT,
        section_type="PUBLIC",
        chart_reference="STEP",
        error=None,
        config=None,
    )

    m = M.MetricInfo(
        column_info=c,
        metric={"data": 1},
        metric_buffers=None,
        metric_summary={"data": 1},
        metric_file_name="1.log",
        metric_step=1,
        metric_epoch=1,
        swanlab_logdir=".",
        swanlab_media_dir=".",
    )
    assert m.column_info.got is None
    assert m.column_info.key == "a/1"
    assert m.column_info.kid == "b"
    assert m.column_info.name == "c"
    assert m.column_info.cls == "SYSTEM"
    assert m.column_info.section_name == "e"
    assert m.column_info.section_sort == 1
    assert m.column_info.chart_type == ChartType.TEXT
    assert m.column_info.chart_reference == "STEP"
    assert m.column_info.section_type == "PUBLIC"
    assert m.column_info.error is None
    assert m.column_info.config is None
    assert m.column_info.key_encode == "a%2F1"
    assert m.column_info.got is None
    assert m.column_info.expected is None
    assert m.column_info.key_encode == "a%2F1"
    assert m.column_info.key == "a/1"
    assert m.column_info.kid == "b"
    assert m.metric == {"data": 1}
    assert m.metric_buffers is None
    assert m.metric_summary == {"data": 1}
    assert m.metric_step == 1
    assert m.metric_epoch == 1
    assert m.swanlab_media_dir == "."
    assert m.metric_file_path == f"./{c.kid}/1.log"
    assert str(m) == "1"
    assert repr(m) == "1"
