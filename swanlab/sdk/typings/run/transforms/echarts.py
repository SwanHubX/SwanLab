"""
@author: cunyue
@file: echarts.py
@time: 2026/4/29
@description: ECharts 模块类型标注
"""

from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from pyecharts.charts.base import Base
    from pyecharts.components.table import Table

    from swanlab.sdk.internal.run.transforms.echarts import ECharts


# 单条数据：ECharts 实例或任意 pyecharts 图表对象
EChartsDataType = Union["Base", "ECharts", "Table"]

# log_echarts 的 data 参数：单条或列表
EChartsDatasType = Union["ECharts", List["ECharts"], "Base", List["Base"], "Table", List["Table"]]
