from ..settings import SwanDataSettings
from ..modules import BaseType, DataType
from ..modules.chart import Chart
from typing import Union, Dict, Protocol
from .utils import create_time, check_key_format, get_a_lock


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, settings: SwanDataSettings) -> None:
        self.settings = settings
        # 当前实验的所有tag数据字段
        self.tags: Dict[str, SwanLabTag] = {}

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
        tag_obj = self.tags.get(tag)
        if tag_obj is None:
            # 添加tag,同时添加图表
            self.tags[tag] = SwanLabTag()
            tag_obj = self.tags[tag]
            tag_obj.create_chart(data)

        # 添加tag信息


class SwanLabTag:
    """
    运行时每个tag的配置，用于记录一些信息
    包括每个tag专属的图表配置
    """

    def __init__(self) -> None:
        self.__sum = 0
        self.__chart = None

    @property
    def sum(self):
        return self.__sum

    @property
    def is_chart_valid(self) -> bool:
        """判断当前tag对应的自动创建图表是否成功
        成功则返回False，失败则返回True
        为True则一切正常
        为False则tag对应的路径不存在
        """
        return False

    def create_chart(self, data: DataType):
        """创建图表，图表可能会创建失败"""
