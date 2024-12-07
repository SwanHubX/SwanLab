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


class HardwareCollector(ABC):
    @abstractmethod
    def collect(self) -> List[HardwareInfo]:
        pass

    def __call__(self):
        try:
            return self.collect()
        except NotImplementedError as n:
            raise n
        except Exception as e:
            swanlog.error("Hardware info collection failed: %s, %s", self.__class__.__name__, str(e))
            return None


# 定义硬件信息执行函数的返回结果
HardwareFuncResult = Tuple[Optional[Any], Optional[HardwareCollector]]
