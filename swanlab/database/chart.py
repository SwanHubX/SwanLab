from .table import ProjectTablePoxy
import ujson
import os
from ..env import swc
from ..utils import create_time, get_a_lock
from typing import List, Union
from ..log import swanlog as swl
from .modules import BaseType


class ChartTable(ProjectTablePoxy):
    """图表管理类，用于管理图表，包括创建图表，删除图表，修改图表配置等操作"""

    default_data = {"_sum": 0, "charts": [], "namespaces": []}

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
        # 当前正在使用的图表
        self.now_charts = []

    def new_chart(
        self,
        chart_id: int,
        namespace: str = "default",
        reference: str = "step",
        chart_type: str = "default",
        config: dict = {},
    ) -> dict:
        """创建一个新chart的配置选项"""
        return {
            "chart_id": chart_id,
            "namespace": namespace,
            "source": [],
            "reference": reference,
            "type": chart_type,
            "config": config,
            "experiment_id": self.experiment_id,
            "update_time": create_time(),
            "create_time": create_time(),
        }

    def add_chart(self, charts, chart):
        """添加图表到配置，同时更新组

        Parameters
        ----------
        data : dict
            配置
        chart : dict
            添加的图表
        """
        namespace: str = chart["namespace"]
        namespaces: list = charts["namespaces"]
        charts["charts"].append(chart)
        # 遍历data["namespaces"]
        ns: dict = None
        exists: bool = False
        for ns in namespaces:
            if ns["namespace"] == namespace:
                exists = True
                break
        # 如果命名空间不存在，添加
        if not exists:
            swl.debug(f"Namespace {namespace} not found, add.")
            ns = {"namespace": namespace, "charts": []}
            if ns["namespace"] == "default":
                swl.debug(f"Namespace {namespace} Add to the beginning")
                namespaces.insert(0, ns)
            else:
                swl.debug(f"Namespace {namespace} Add to the end.")
                namespaces.append(ns)
        # 添加当前的chart_id到结尾
        ns["charts"].append(chart["chart_id"])
        swl.debug(f"Chart {chart['chart_id']} add, now charts: " + str(ns["charts"]))

    def add(self, tag: str, data: Union[float, int, BaseType]):
        """添加一张图表

        Parameters
        ----------
        tag : str
            图表标签
        data: Union[float, int, BaseType]
            数据源，可以是一个数字，也可以是一个swanlab.BaseType的子类
        """
        with get_a_lock(self.path) as f:
            charts = ujson.load(f)
            # 记录图表数量，+1
            charts["_sum"] += 1
            chart = self.new_chart(charts["_sum"])
            # 如果data是BaseType类型，解构，并且修改一些必要参数
            if issubclass(type(data), BaseType):
                chart["namespace"], (chart["type"], data_types), chart["reference"], chart["config"] = data.__next__()
            else:
                data_types = [float, int]
            # 如果data不是data_types中的类型，尝试转换为这两个类型中的一个（优先转换为float）
            if not isinstance(data, tuple(data_types)):
                try:
                    data = self.try_convert(data, data_types)
                except:
                    # 此时代表数据异常，拿到data的__class__.__name__，记录到chart.error中并保存
                    class_name = data.__class__.__name__
                    excepted = [i.__name__ for i in data_types]
                    swl.error(f"Data type error, tag: {tag}, data type: {class_name}, excepted: {excepted}")
                    chart["error"] = {"data_class": class_name, "excepted": excepted}
            chart["source"].append(tag)
            # 添加图表
            self.add_chart(charts, chart)
            f.truncate(0)
            f.seek(0)
            ujson.dump(charts, f, ensure_ascii=False)
            # 记录当前chart期望的数据类型，这个字段不会被保存到文件中
            chart["excepted"] = data_types
            self.now_charts.append(chart)
            f.close()

    def is_chart_error(self, tag):
        """遍历所有自动创建的chart，检查对应的tag的chart是否存在error字段

        Parameters
        ----------
        tag : str
            tag名称
        """
        for chart in self.now_charts:
            if tag in chart["source"] and "error" in chart:
                return chart["error"]
        return False

    def try_convert(self, value, data_types):
        # 如果当前data已经是data_types中的类型，直接返回
        if isinstance(value, tuple(data_types)):
            return value
        for data_type in data_types:
            try:
                converted_value = data_type(value)
                return converted_value
            except (ValueError, TypeError):
                # 如果转换失败，继续尝试下一种类型
                continue

        # 如果所有类型都尝试过仍然失败，则抛出异常
        raise ValueError(f"Unable to convert {value} to any of the specified types.")

    def try_convert_data(self, data, tag):
        """尝试将data转换为data_types中的类型，如果转换失败，返回None
        如果所有类型都尝试过仍然失败，则抛出异常
        调用这个函数代表用户的图表已经创建并且传入了与之前不同的数据类型
        """
        # 首先寻找tag对应的chart
        for chart in self.now_charts:
            if tag in chart["source"]:
                data_types = chart["excepted"]
                break
        # 如果当前data已经是data_types中的类型，直接返回
        if isinstance(data, tuple(data_types)):
            return data
        try:
            return self.try_convert(data, data_types)
        except:
            raise ValueError(
                f"Unable to convert {data} to any of the specified types. This Error means that you have changed the data type of the chart, please check the data type of the chart."
            )
