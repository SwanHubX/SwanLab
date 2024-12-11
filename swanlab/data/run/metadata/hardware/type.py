"""
@author: cunyue
@file: type.py
@time: 2024/12/3 20:27
@description: 硬件信息采集类型定义
"""

from abc import ABC, abstractmethod
from typing import TypedDict, Tuple, Optional, Any, List, Union

from swankit.callback.models import ColumnConfig

from swanlab.log import swanlog


# 定义硬件信息类型
class HardwareInfo(TypedDict):
    # 硬件信息键
    key: str
    # 硬件信息值
    value: Union[str, int, float]
    # 硬件信息名称
    name: str
    # 相关配置
    config: Optional[ColumnConfig]


HardwareInfoList = List[HardwareInfo]


class CollectGuard:
    """
    采集守卫，在采集任务执行前后可以选择执行一些操作
    这些操作可能比较耗费资源，因此需要选择性执行，而选择的依据是是执行的次数
    这与任务的采集间隔有关，即此类适配了 MonitorCron 的采集间隔
    """

    def __init__(self):
        self.collect_num = 0

    def before_collect(self):
        try:
            # count=0，执行一次
            if self.collect_num == 0:
                return self.before_collect_impl()
            # 60>count>0, 不执行
            elif self.collect_num < 60:
                return
            else:
                return self.before_collect_impl()
        finally:
            self.collect_num += 1

    def after_collect(self):
        # count=1, 执行一次
        if self.collect_num == 1:
            return self.after_collect_impl()
        # 61>count>1, 不执行
        elif self.collect_num < 61:
            return
        else:
            return self.after_collect_impl()

    def before_collect_impl(self):
        pass

    def after_collect_impl(self):
        pass


class HardwareCollector(CollectGuard, ABC):

    @abstractmethod
    def collect(self) -> HardwareInfoList:
        pass

    def __call__(self) -> Optional[HardwareInfoList]:
        try:
            self.before_collect()
            self.collect()
            self.after_collect()
            return
        except NotImplementedError as n:
            raise n
        except Exception as e:
            swanlog.error(f"Hardware info collection failed: {self.__class__.__name__}, {str(e)}")
            return None


# 定义硬件信息执行函数的返回结果
HardwareFuncResult = Tuple[Optional[Any], Optional[HardwareCollector]]
