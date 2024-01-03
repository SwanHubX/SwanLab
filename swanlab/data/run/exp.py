from ..settings import SwanDataSettings
from ..modules import BaseType, DataType
from ...log import swanlog
from typing import Dict
from .utils import create_time, check_key_format, get_a_lock
from urllib.parse import quote
import ujson
import os
import math


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, settings: SwanDataSettings, id: int) -> None:
        self.settings = settings
        if not os.path.exists(self.settings.log_dir):
            os.mkdir(self.settings.log_dir)
        # 当前实验的所有tag数据字段
        self.tags: Dict[str, SwanLabTag] = {}
        self.id = id
        if not os.path.exists(self.settings.chart_path):
            with open(self.settings.chart_path, "w", encoding="utf-8") as f:
                f.write(ujson.dumps(self.__new_charts()))

    def add(self, tag: str, data: DataType, step: int = None):
        """记录一条新的tag数据

        Parameters
        ----------
        tag : str
            tag名称，需要检查数据类型
        data : Union[int, float, BaseType]
            tag数据，可以是浮点型，也可以是SwanLab定义的数据类型，具体与添加的图表类型有关系
            如果data是int或者float，添加图表时自动添加默认折线图，如果是BaseType，添加图表时选择对应的图表
            需要检查数据类型

        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        check_key_format(tag)
        if isinstance(data, BaseType):
            # 注入一些内容
            data.settings = self.settings
        tag_obj = self.tags.get(tag)
        if tag_obj is None:
            # 添加tag,同时添加图表
            tag_obj = SwanLabTag(self.id, tag, self.settings.log_dir, self.settings.chart_path)
            self.tags[tag] = tag_obj
            """
            由于添加图表的同时会尝试转换data的类型，但这在注入step之前
            所以此处需要手动注入依赖
            关于step，如果step格式不正确，会在add的时候被拦截，此处如果格式不正确直接设置为1即可
            """
            if isinstance(data, BaseType):
                data.step = 1 if step is None or not isinstance(step, int) else step
                data.tag = tag_obj.tag
            tag_obj.create_chart(tag, data)
        if not tag_obj.is_chart_valid:
            return swanlog.warning(f"Chart {tag} has been marked as error, ignored.")
        # 添加tag信息
        tag_obj.add(data, step)

    def __new_charts(self):
        return {"_sum": 0, "charts": [], "namespaces": []}


class SwanLabTag:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

    def __init__(self, experiment_id, tag: str, log_dir: str, chart_path: str) -> None:
        self.experiment_id = experiment_id
        self.tag = tag
        self.__steps = set()
        # 当前tag的图表配置
        self.__chart = None
        self.__log_dir = log_dir
        self.__chart_path = chart_path
        # 默认数据类型
        self.data_types = [float, int]
        # summary 数据概要总结
        self._summary = {}
        # 当前tag的存储数据
        self.__data = self.__new_tags()

    @property
    def sum(self):
        return len(self.__steps)

    @property
    def is_chart_valid(self) -> bool:
        """判断当前tag对应的自动创建图表是否成功
        成功则返回False，失败则返回True
        为True则一切正常
        为False则tag对应的路径不存在
        """
        return self.__chart.get("error") is None

    def add(self, data: DataType, step: int = None):
        """添加一个数据，在内部完成数据类型转换
        如果转换失败，打印警告并退出
        并且添加数据，当前的数据保存是直接保存，后面会改成缓存形式

        Parameters
        ----------
        data : DataType
            待添加的数据
        """
        # 如果step不是None也不是int，设置为None且打印警告
        if step is not None and not isinstance(step, int):
            swanlog.warning(f"Step {step} is not int, SwanLab will set it automatically.")
            step = None
        # 转换step，如果step为None，则改为len(self.__steps)+1
        if step is None:
            step = len(self.__steps) + 1
        # 如果step已经存在，打印警告并退出
        if step in self.__steps:
            return swanlog.warning(f"Step {step} on tag {self.tag} already exists, ignored.")
        # 添加数据，首先比较数据，如果数据比之前的数据大，则更新最大值，否则不更新
        data = self.try_convert_data(data, step)
        # 如果数据比之前的数据小，则更新最小值，否则不更新
        self._summary["max"] = data if self._summary.get("max") is None else max(self._summary["max"], data)
        self._summary["min"] = data if self._summary.get("min") is None else min(self._summary["min"], data)
        self._summary["num"] = self._summary.get("num", 0) + 1
        self.__steps.add(step)
        swanlog.debug(f"Add data, tag: {self.tag}, step: {step}, data: {data}")
        # ---------------------------------- 保存数据 ----------------------------------
        self.__add_data(data, step)

    def __add_data(self, data, step: int):
        """添加数据到data中"""
        if len(self.__data["data"]) >= self.__slice_size:
            # 如果当前数据已经达到了__slice_size，重新创建一个新的data
            self.__data = self.__new_tags()
        # 添加数据
        self.__data["data"].append(self.__new_tag(step, data))
        # 优化文件分片，每__slice_size个tag数据保存为一个文件，通过sum来判断
        sum = len(self.__steps)
        mu = math.ceil(sum / self.__slice_size)
        # 存储路径
        file_path = os.path.join(self.save_path, str(mu * self.__slice_size) + ".json")
        # 方便一些，直接使用w+模式覆盖写入
        with get_a_lock(file_path, mode="w+") as f:
            ujson.dump(self.__data, f, ensure_ascii=False)
        # 更新实验信息总结
        with get_a_lock(os.path.join(self.save_path, "_summary.json"), "w+") as f:
            ujson.dump(self._summary, f, ensure_ascii=False)

    @property
    def save_path(self):
        """获取当前tag的保存路径

        Returns
        -------
        str
            保存路径
        """
        path = os.path.join(self.__log_dir, quote(self.tag, safe=""))
        if not os.path.exists(path):
            os.mkdir(path)
        return path

    def create_chart(self, tag: str, data: DataType):
        """创建图表，图表可能会创建失败"""
        with get_a_lock(self.__chart_path) as f:
            charts = ujson.load(f)
            # 记录图表数量，+1
            charts["_sum"] += 1
            chart = self.__new_chart(charts["_sum"], self.experiment_id)
            # 如果data是BaseType类型，解构，并且修改一些必要参数
            if issubclass(type(data), BaseType):
                chart["namespace"], types, chart["reference"], chart["config"] = data.__next__()
                chart["type"], self.data_types = types
            # 如果data不是DataType中的类型，尝试转换为这两个类型中的一个（优先转换为float）
            if not isinstance(data, tuple(self.data_types)):
                try:
                    data = self.try_convert(data)
                except:
                    # 此时代表数据异常，拿到data的__class__.__name__，记录到chart.error中并保存
                    class_name = data.__class__.__name__
                    excepted = [i.__name__ for i in self.data_types]
                    swanlog.error(f"Data type error, tag: {tag}, data type: {class_name}, excepted: {excepted}")
                    chart["error"] = {"data_class": class_name, "excepted": excepted}
            chart["source"].append(tag)
            # 添加图表
            self.__add_chart(charts, chart)
            f.truncate(0)
            f.seek(0)
            ujson.dump(charts, f, ensure_ascii=False)

    def __add_chart(self, charts, chart):
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
            swanlog.debug(f"Namespace {namespace} not found, add.")
            ns = {"namespace": namespace, "charts": []}
            if ns["namespace"] == "default":
                swanlog.debug(f"Namespace {namespace} Add to the beginning")
                namespaces.insert(0, ns)
            else:
                swanlog.debug(f"Namespace {namespace} Add to the end.")
                namespaces.append(ns)
        # 添加当前的chart_id到结尾
        ns["charts"].append(chart["chart_id"])
        swanlog.debug(f"Chart {chart['chart_id']} add, now charts: " + str(ns["charts"]))
        self.__chart = chart

    def __new_chart(self, chart_id, experiment_id):
        return {
            "chart_id": chart_id,
            "namespace": "default",
            "source": [],
            "reference": "step",
            "type": "default",
            "config": {},
            "experiment_id": experiment_id,
            "update_time": create_time(),
            "create_time": create_time(),
        }

    def __new_tag(self, index, data) -> dict:
        """创建一个新的data数据，实际上是一个字典，包含一些默认信息"""
        return {
            "index": str(index),
            "data": data,
            "create_time": create_time(),
        }

    def __new_tags(self) -> dict:
        """创建一个新的tag data数据集合

        Returns
        -------
        dict
            返回一个新的data数据集合
        """
        time = create_time()
        return {
            "create_time": time,
            "update_time": time,
            "data": [],
        }

    def try_convert(self, value: DataType):
        # 如果当前data已经是data_types中的类型，直接返回
        if isinstance(value, tuple(self.data_types)):
            return value
        if isinstance(value, BaseType):
            value = value.convert
        for data_type in self.data_types:
            try:
                converted_value = data_type(value)
                return converted_value
            except (ValueError, TypeError):
                # 如果转换失败，继续尝试下一种类型
                continue
        # 如果所有类型都尝试过仍然失败，则抛出异常
        raise ValueError(f"Unable to convert {value} to any of the specified types.")

    def try_convert_data(self, data, step):
        """尝试将data转换为data_types中的类型，如果转换失败，返回None
        如果所有类型都尝试过仍然失败，则抛出异常
        调用这个函数代表用户的图表已经创建并且传入了与之前不同的数据类型
        """
        try:
            if isinstance(data, BaseType):
                # 注入内容
                if data.step is None:
                    data.step = step
                if data.tag is None:
                    data.tag = self.tag
            return self.try_convert(data)
        except ValueError:
            raise ValueError(
                f"Unable to convert {data} to any of the specified types. This Error means that you have changed the data type of the chart, please check the data type of the chart."
            )
