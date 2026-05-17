"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13

@description: SwanLab Core Python 版本，封装SwanLab云端版核心业务，包括：
1. 提供http客户端，用于与SwanLab云端API进行交互。
2. 提供rpc封装函数，以rpc方式调用SwanLab云端API。
3. 提供上传线程，在另一个线程执行上传任务。
4. 存储指标上下文，便于分布式训练时的指标同步。
...

实现 CoreProtocol，当前为纯 Python 实现。
未来由 swanlab-core（Go 二进制）替代时，此模块整体被替换，
BackgroundConsumer 等调用方无需修改。


Core 同时需要根据不同模式处理不同的业务，这是设计模式决定的
值得说明的是，在当前的上层设计中，upsert 方法在 disabled 模式下永远不会触发，但是考虑到设计完整性，我们增加了相关业务逻辑判断
"""

from .core import CorePython
from .sync import CoreSyncPython

__all__ = ["CorePython", "CoreSyncPython"]
