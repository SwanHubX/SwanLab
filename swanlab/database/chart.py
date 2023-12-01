from .table import ProjectTablePoxy
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

    def __init__(self, base_path: str):
        """初始化图表管理类"""
        path = os.path.join(base_path, "charts.json")
        # 判断path是否存在，如果存在，则加载数据，否则创建
        if os.path.exists(path):
            with open(path, "r") as f:
                data = ujson.load(f)
        else:
            data = {
                "charts": [],
            }
        # 保存表单信息，创建时间等内容在父类中完成
        super().__init__(data, path)
        # 保存表单信息
        self.save()
