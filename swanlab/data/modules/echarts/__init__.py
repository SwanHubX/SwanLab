import json
from swankit.core import MediaBuffer, DataSuite as D, MediaType
from pyecharts.charts.base import Base

class Echarts(MediaType):
    def __init__(self, chart: Base):
        self._chart = chart
     # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        options = self._chart.dump_options()
        filename = f"echart_{id(self)}.json"
        buffer = json.dumps(options).encode("utf-8")
        return filename, MediaBuffer()

    
    def get_chart(self):
        return self.Chart.ECHARTS

    
    def get_section(self):
        return "ECharts"