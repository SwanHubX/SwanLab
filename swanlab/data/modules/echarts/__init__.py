from swankit.core import MediaBuffer, DataSuite as D, MediaType
from pyecharts.charts.base import Base

class Echarts(MediaType):
    def __init__(self, chart: Base):
        self._chart = chart
        self.buffer = MediaBuffer()
        self.buffer.write(self._chart.dump_options().encode('utf-8'))  # 写入二进制数据
        self.buffer.seek(0)  # 重置指针到开头，以便后续读取
     # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        hash_name = D.get_hash_by_bytes(self._chart.dump_options().encode('utf-8'))[:16]
        filename = f"echart-step{self.step}-{hash_name}.json"
        return filename,self.buffer

    def get_chart(self):
        return self.Chart.ECHARTS

    def get_section(self):
        return "ECharts"