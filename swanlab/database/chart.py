from .table import ProjectTablePoxy
import ujson
import os
from ..utils import create_time


class ChartTable(ProjectTablePoxy):
    """图表管理类，用于管理图表，包括创建图表，删除图表，修改图表配置等操作"""

    default_data = {"_sum": 0, "charts": []}

    def __init__(self, base_path: str, experiment_id: int):
        """初始化图表管理类"""
        # 判断path是否存在，如果存在，则加载数据，否则创建
        self.experiment_id = experiment_id
        self.path = os.path.join(base_path, "charts.json")
        if os.path.exists(self.path):
            with open(self.path, "r") as f:
                data = ujson.load(f)
        else:
            data = self.default_data
        # 保存表单信息，创建时间等内容在父类中完成
        super().__init__(data, self.path)
        # 保存表单信息
        self.save_no_lock()

    def new_chart(
        self, chart_id: int, namespace: str = "default", reference: str = "step", chart_type: str = "default"
    ) -> dict:
        """创建一个新chart的配置选项"""
        return {
            "chart_id": chart_id,
            "namespace": namespace,
            "source": [],
            "reference": reference,
            "type": chart_type,
            "config": {},
            "experiment_id": self.experiment_id,
            "update_time": create_time(),
            "create_time": create_time(),
        }

    def add(self, tag: str, namespace: str = None, reference: str = None, chart_type: str = None, config: dict = None):
        """添加一个图标"""
        data = {}
        with open(self.path, "r") as f:
            data = ujson.load(f)
            # 记录数据
            data["_sum"] += 1
            chart = self.new_chart(chart_id=data["_sum"])
            chart["source"].append(tag)
            data["charts"].append(chart)
        with open(self.path, "w") as f:
            ujson.dump(data, f)
        pass
