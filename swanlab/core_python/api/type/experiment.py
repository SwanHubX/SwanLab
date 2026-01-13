"""
@author: Zhou QiYang
@file: type.py
@time: 2026/1/10 22:12
@description: 实验相关后端API接口类型
"""

from typing import TypedDict, List, Dict, Optional, Literal

ColumnType = Literal['STABLE', 'SCALAR', 'CONFIG']  # 列类型
StateType = Literal['FINISHED', 'CRASHED', 'ABORTED', 'RUNNING']  # 实验状态


# ------------------------------------- 通用类型 -------------------------------------
class UserType(TypedDict):
    username: str  # 用户名
    name: str  # 用户显示名称


class RunType(TypedDict):
    cuid: str  # 实验CUID, 唯一标识符
    name: str  # 实验名称
    createdAt: str  # 创建时间, e.g., '2024-11-23T12:28:04.286Z'
    description: str  # 实验描述
    labels: List[Dict[str, str]]  # 实验标签列表
    profile: Dict[str, Dict[str, object]]  # 实验配置和摘要信息，包含 'config' 和 'scalar'
    state: StateType  # 实验状态
    cluster: str  # 实验组
    job: str  # 任务类型
    runtime: str  # 运行时间
    user: UserType  # 实验所属用户
    rootExpId: Optional[str]  # 祖宗实验对应的实验 cuid，如果为克隆实验则必传
    rootProId: Optional[str]  # 祖宗实验对应的项目 cuid，如果为克隆实验则必传
