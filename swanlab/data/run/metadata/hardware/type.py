"""
@author: cunyue
@file: type.py
@time: 2024/12/3 20:27
@description: 硬件信息采集类型定义
"""

from typing import TypedDict, Callable, Tuple, Optional, Any, List, Union


# 定义硬件信息类型
class HardwareInfo(TypedDict):
    # 硬件信息键
    key: str
    # 硬件信息值
    value: Union[str, int, float]
    # 硬件信息名称
    name: str


# 定义硬件信息采集函数类型
HardwareMonitorFunc = Callable[[], Optional[HardwareInfo]]

# 定义硬件信息执行函数的返回结果
HardwareFuncResult = Tuple[Optional[Any], List[HardwareMonitorFunc]]
