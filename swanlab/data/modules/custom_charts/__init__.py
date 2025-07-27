"""
@author: ComPleHN
@file: __init__.py
@time: 2025/5/19 14:01
@desc: 自定义图表，目前集成了 echarts
"""

from typing import Union

import pyecharts
from pyecharts.charts.base import Base

from swanlab.toolkit import MediaBuffer, DataSuite as D, MediaType
from . import echarts
from .table import Table
from .metrics import confusion_matrix, pr_curve, roc_curve

PyEchartsBase = pyecharts.charts.base.Base
"""
pyecharts.charts.base.Base
"""
PyEchartsTable = Table
"""
custom Table inherited from pyecharts.components.table.Table
"""

__all__ = ["echarts", 'Echarts', 'PyEchartsTable', 'PyEchartsBase', "roc_curve", "pr_curve", "confusion_matrix"]


class Echarts(MediaType):
    def __init__(self, chart: Union[PyEchartsBase, PyEchartsTable]):
        super().__init__()
        self._chart = chart
        self.buffer = MediaBuffer()

    # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        byte_string = self._chart.dump_options().encode('utf-8')
        hash_name = D.get_hash_by_bytes(byte_string)[:16]

        # 写入buffer
        self.buffer.write(byte_string)  # 写入二进制数据
        self.buffer.seek(0)  # 重置指针到开头，以便后续读取

        filename = f"echarts-step{self.step}-{hash_name}.json"
        return filename, self.buffer

    def get_chart(self):
        return self.Chart.ECHARTS

    def get_section(self):
        return "ECharts"
