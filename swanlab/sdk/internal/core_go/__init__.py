"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: SwanLab Go core 的 Python RPC 客户端。

CoreEnum.CORE（"Core"）对应本包 CoreGo；CoreEnum.CORE_PYTHON 对应 CorePython。
公共出口只做 re-export，不承载逻辑。
"""

from swanlab.sdk.internal.core_go.core import CoreGo
from swanlab.sdk.internal.core_go.errors import CoreGoError, ReasonCode, Stage, reason_for_grpc

__all__ = ["CoreGo", "CoreGoError", "ReasonCode", "Stage", "reason_for_grpc"]
