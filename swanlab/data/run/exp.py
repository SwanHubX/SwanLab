from ..settings import SwanDataSettings
from ..modules import BaseType, DataType
from ...log import swanlog
from typing import Dict
from .utils import create_time, check_tag_format, get_a_lock
from urllib.parse import quote
import ujson
import os
import math
from datetime import datetime
from .db import (
    Tag,
    Namespace,
    Chart,
    Source,
    Experiment,
    Display,
    ExistedError,
)


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, settings: SwanDataSettings, id: int, exp: Experiment) -> None:
        """初始化实验

        Parameters
        ----------
        settings : SwanDataSettings
            全局运行时配置
        id : int
            实验id
        exp : Experiment
            数据库实例，代表当前实验行
        """
        self.settings = settings
        if not os.path.exists(self.settings.log_dir):
            os.mkdir(self.settings.log_dir)
        # 当前实验的所有tag数据字段
        self.tags: Dict[str, SwanLabTag] = {}
        self.id = id
        """此实验对应的id，实际上就是db.id"""
        self.db = exp
        """此实验对应的数据库实例"""

    def add(self, key: str, data: DataType, step: int = None):
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
        tag = check_tag_format(key, auto_cut=True)
        if key != tag:
            # 超过255字符，截断
            swanlog.warning(f"Tag {key} is too long, cut to 255 characters.")

        # 如果是swanlab自定义的数据类型，注入settings，方便在类内部使用
        if isinstance(data, BaseType):
            data.settings = self.settings
        # 判断tag是否存在，如果不存在则创建tag
        tag_obj = self.tags.get(tag)
        """
        数据库创建字段将在chart创建完成后进行
        无论格式是否正确，都会创建chart，但是如果格式不正确，不会写入日志
        """
        if tag_obj is None:
            # 将此tag对象添加到实验列表中
            tag_obj = SwanLabTag(self.id, tag, self.settings.log_dir)
            self.tags[tag] = tag_obj
            """
            由于添加图表的同时会尝试转换data的类型，但这在注入step之前
            所以此处需要手动注入依赖
            关于step，如果step格式不正确，会在tag_obj.add的时候被拦截，此处如果格式不正确直接设置为0即可
            """
            if isinstance(data, BaseType):
                data.step = 0 if step is None or not isinstance(step, int) else step
                data.tag = tag_obj.tag
            tag_obj.create_chart(tag, data)
        # 检查tag创建时图表是否创建成功，如果失败则也没有写入数据的必要了，直接退出
        if not tag_obj.is_chart_valid:
            return swanlog.warning(f"Chart {tag} has been marked as error, ignored.")
        # 添加tag信息
        step = tag_obj.add(data, step)
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        swanlog.log(f"{timestamp}  step:{step}  tag:{tag}  value:{str(data)}")


class SwanLabTag:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

    def __init__(self, experiment_id, tag: str, log_dir: str) -> None:
        self.experiment_id = experiment_id
        self.tag = tag
        self.__steps = set()
        """此tag已经包含的steps步骤
        """
        self.__chart = None
        """当前tag的图表配置"""
        self.__log_dir = log_dir
        """存储文件夹路径"""
        self.data_types = [float, int]
        """默认数据类型，如果tag数据为BaseType的子类，则使用其规定的数据类型"""
        self._summary = {}
        """数据概要总结"""
        self.__data = self.__new_tags()
        """当前tag的数据"""
        self.__namespace = None
        """当前tag自动生成图表的命名空间"""
        self.__error = None
        """此tag在自动生成chart的时候的错误信息"""

    @property
    def sum(self):
        """当前tag的数据总数"""
        return len(self.__steps)

    @property
    def is_chart_valid(self) -> bool:
        """判断当前tag对应的自动创建图表是否成功
        成功则返回False，失败则返回True
        为True则一切正常
        为False则tag对应的路径不存在
        """
        return self.__error is None

    def __is_nan(self, data):
        """判断data是否为nan"""
        return isinstance(data, (int, float)) and math.isnan(data)

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
        # 转换step，如果step为None，则改为合适的长度，目前step从0开始
        if step is None:
            step = len(self.__steps)
        # 如果step已经存在，打印警告并退出
        if step in self.__steps:
            return swanlog.warning(f"Step {step} on tag {self.tag} already exists, ignored.")
        # 添加数据，首先比较数据，如果数据比之前的数据大，则更新最大值，否则不更新
        """
        python环境下data可能是列表、字符串、整型、浮点型等
        因此下面的summary还需要在未来做好兼容
        对于整型和浮点型，还存在极大和极小值的问题
        目前的策略是让python解释器自己处理，在前端完成数据的格式化展示
        """
        try:
            data = self.try_convert_after_add_chart(data, step)
        except ValueError:
            return swanlog.warning(
                f"Data {data} on tag {self.tag} cannot be converted, SwanLab will ignore it, but the chart still exists."
            )
        is_nan = self.__is_nan(data)
        if not is_nan:
            # 如果数据比之前的数据小，则更新最小值，否则不更新
            self._summary["max"] = data if self._summary.get("max") is None else max(self._summary["max"], data)
            self._summary["min"] = data if self._summary.get("min") is None else min(self._summary["min"], data)
        self._summary["num"] = self._summary.get("num", 0) + 1
        self.__steps.add(step)
        swanlog.debug(f"Add data, tag: {self.tag}, step: {step}, data: {data}")
        # ---------------------------------- 保存数据 ----------------------------------
        """添加数据到data中"""
        if len(self.__data["data"]) >= self.__slice_size:
            # 如果当前数据已经达到了__slice_size，重新创建一个新的data
            self.__data = self.__new_tags()
        # 添加数据
        data = data if not is_nan else "NaN"
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

        return step

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
        """在第一次添加tag的时候，自动创建图表和namespaces，同时写入tag和数据库信息，将创建的信息保存到数据库中
        此方法只能执行一次
        具体步骤是：
        1. 创建tag字段
        2. 如果不是baseType类型，

        """
        if self.__namespace is not None or self.__chart is not None:
            raise ValueError(f"Chart {tag} has been created, cannot create again.")

        # 如果是非BaseType类型，写入默认命名空间，否则写入BaseType指定的命名空间
        if not isinstance(data, BaseType):
            # data_type不变
            namespace, chart_type, reference, config = "default", "default", "step", None
            sort = 0
        else:
            # 解构数据
            namespace, types, reference, config = data.__next__()
            chart_type, self.data_types = types
            sort = None
        # 创建chart
        self.__chart: Chart = Chart.create(
            tag,
            experiment_id=self.experiment_id,
            type=chart_type,
            reference=reference,
            config=config,
        )
        # 创建命名空间，如果命名空间已经存在，会抛出ExistedError异常，捕获不处理即可
        # 需要指定sort，default命名空间的sort为0，其他命名空间的sort为None，表示默认添加到最后
        try:
            self.__namespace = Namespace.create(name=namespace, experiment_id=self.experiment_id, sort=sort)
            swanlog.debug(f"Namespace {namespace} created, id: {self.__namespace.id}")
        except ExistedError:
            self.__namespace: Namespace = Namespace.get(name=namespace, experiment_id=self.experiment_id)
            swanlog.debug(f"Namespace {namespace} exists, id: {self.__namespace.id}")
        # 创建display，这个必然是成功的，因为display是唯一的，直接添加到最后一条即可
        Display.create(chart_id=self.__chart.id, namespace_id=self.__namespace.id)
        """
        接下来判断tag格式的正确性，判断完毕后往source中添加一条tag记录
        在此函数中，只判断tag的格式是否正确，不记录数据
        在逻辑上只有第一次会检查tag的正确性，也就是说前端的error错误只有在第一次添加tag的时候才有可能出现
        如果第一次添加成功，后续出现错误，只会在添加的时候warning一下然后丢弃这个错误
        如果第一次添加失败，后续都不会再添加此数据，因此不会出现错误
        """
        # 如果data不是期望的data_types中的类型，尝试转换为这两个类型中的一个（优先转换为第一个）
        # 如果data是BaseType类型，会在try_convert中完成转换，此处不需要管
        error = None
        if not isinstance(data, tuple(self.data_types)):
            try:
                data = self.try_convert(data)
            except:
                # 此时代表数据异常，拿到data的__class__.__name__，生成error并保存
                class_name = data.__class__.__name__
                excepted = [i.__name__ for i in self.data_types]
                swanlog.error(f"Data type error, tag: {tag}, data type: {class_name}, excepted: {excepted}")
                error = {"data_class": class_name, "excepted": excepted}
        if self.__is_nan(data):
            """如果data是nan，生成error并保存"""
            error = {"data_class": "NaN", "excepted": [i.__name__ for i in self.data_types]}
        # 添加一条tag记录
        tag: Tag = Tag.create(
            experiment_id=self.experiment_id,
            name=tag,
            # TODO 类型暂时为default
            type="default",
        )
        # 添加一条source记录
        Source.create(tag_id=tag.id, chart_id=self.__chart.id, error=error)
        self.__error = error

    def __new_tag(self, index, data) -> dict:
        """创建一个新的data数据，实际上是一个字典，包含一些默认信息"""
        return {
            "index": str(index),
            "data": data,
            "create_time": create_time(),
        }

    @staticmethod
    def __new_tags() -> dict:
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

    def try_convert_after_add_chart(self, data, step):
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
