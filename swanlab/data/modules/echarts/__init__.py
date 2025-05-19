import json
from swankit.core import MediaBuffer, DataSuite as D, MediaType
from pyecharts.charts.base import Base

class Echarts(MediaType):
    def __init__(self, chart: Base):
        self._chart = chart
     # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        # hash_name = D.get_hash_by_ndarray(self._chart.dump_options())[:16]
        hash_name = '123'
        filename = f"echart-{hash_name}.json"
        return filename, MediaBuffer()

    def get_chart(self):
        return self.Chart.ECHARTS

    def get_section(self):
        return "ECharts"