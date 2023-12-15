from .table import ProjectTablePoxy
import ujson
import os
from ..env import swc
from ..utils import create_time, get_a_lock
from typing import List, Union


class ChartTable(ProjectTablePoxy):
    """图表管理类，用于管理图表，包括创建图表，删除图表，修改图表配置等操作"""

    default_data = {"_sum": 0, "charts": []}

    def __init__(self, experiment_id: int):
        """初始化图表管理类"""
        # 判断path是否存在，如果存在，则加载数据，否则创建
        self.experiment_id = experiment_id
        # 文件保存路径
        self.path = swc.chart
        if os.path.exists(self.path):
            with open(self.path, "r", encoding="utf-8") as f:
                data = ujson.load(f)
        else:
            data = self.default_data
        # 保存表单信息，创建时间等内容在父类中完成
        super().__init__(data, self.path)
        # 保存表单信息
        self.save_no_lock()

    def new_chart(self, chart_id: int, namespace: str, reference: str, chart_type: str) -> dict:
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

    def add(
        self,
        tag: Union[str, List[str]],
        namespace: str = "default",
        reference: str = "step",
        chart_type: str = "default",
    ):
        """添加一张图表

        Parameters
        ----------
        tag : Union[str, List[str]]
            图表标签，可以是一个标签，也可以是多个标签，但必须是字符串或者字符串列表
            单标签代表创建的图表只包含一个数据源，多标签代表创建的图表包含多个数据源
        namespace : str, optional
            命名空间，用于区分不同的图表在前端的显示位置
        reference : str, optional
            参考系，用于区分不同的图表表格的组织方式
        chart_type : str, optional
            图表类型，用于区分不同的图表的显示方式，如折线图，柱状图等
        """
        f = get_a_lock(self.path)
        data = ujson.load(f)
        # 记录数据
        data["_sum"] += 1
        chart = self.new_chart(data["_sum"], namespace, reference, chart_type)
        chart["source"].append(tag)
        data["charts"].append(chart)
        f.truncate(0)
        f.seek(0)
        ujson.dump(data, f, ensure_ascii=False)
        f.close()
