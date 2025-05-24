"""
@author: ComPleHN
@file: __init__.py
@time: 2025/5/19 14:01
@desc: 集成 pyecharts
"""

import pyecharts
from pyecharts.charts.base import Base
from swankit.core import MediaBuffer, DataSuite as D, MediaType

echarts = pyecharts.charts

PyEchartsBase = pyecharts.base.Base  # noqa: 无法解析.base文件的导入，但实际上可以导入
"""
pyecharts.base.Base的别名
"""

__all__ = ["echarts", 'Echarts', 'PyEchartsBase']


class Echarts(MediaType):
    def __init__(self, chart: Base):
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
