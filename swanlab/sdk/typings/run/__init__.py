"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时类型
"""

from typing import Literal

ResumeType = Literal["must", "allow", "never"]
"""
Run 恢复策略
"""

ParallelType = Literal["none", "shared"]
"""
Run 并行策略
"""

ModeType = Literal["disabled", "online", "local", "offline"]
"""
Run 运行策略
"""

OnlineErrorPolicy = Literal["disabled", "offline", "raise"]
"""
在线模式连接失败策略

仅连接类 / 服务端瞬时错误参与降级；认证与配置错误始终抛错。

* `disabled`: 服务不可用时降级为 disabled，训练继续但不记录数据。
* `offline`: 服务不可用时降级为 offline，继续本地记录。
* `raise`: 保持严格失败，连接异常时直接抛错。
"""

FinishType = Literal["success", "crashed", "aborted"]
"""
Run 结束策略
"""


AsyncLogType = Literal["asyncio", "threading", "spawn", "fork"]
"""
异步日志上报类型
"""

SaveType = Literal["live", "now", "end"]
"""
文件保存策略
"""
