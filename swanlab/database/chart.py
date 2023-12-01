from .table import ProjectTablePoxy
from ..env import SWANLAB_LOGS_FOLDER
import ujson
import os

DEFAULT_CHART = {
    "chart_id": 0,
    "tag": "default",
    "source": [],
    "type": "default",
    "config": {},
}


class ChartTable(ProjectTablePoxy):
    """图表管理类，用于管理图表，包括创建图表，删除图表，修改图表配置等操作"""

    path = os.path.join(SWANLAB_LOGS_FOLDER, "chart.json")

    def __init__(self):
        """初始化图表管理类"""
        # 判断path是否存在，如果存在，则加载数据，否则创建
        if os.path.exists(self.path):
            with open(self.path, "r") as f:
                data = ujson.load(f)
        else:
            data = {
                "charts": [],
            }
        # 保存表单信息
        super().__init__(data, self.path)
        # 保存表单信息
        self.save()
