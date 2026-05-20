"""
@author: cunyue
@file: __init__.py
@time: 2026/3/31 01:44
@description: 硬件提供商模块
由于硬件同时包含多维度信息，也同时涉及到信息采集和监控，不同厂家的硬件信息采集和监控方式也不一样，所以单独成一个模块
模块下每个子模块负责一个厂商的硬件信息采集和监控，模块内提供统一接口供上层调用

设计取舍 —— collect() 的两种实现模式：

CollectorProtocol 提供了基于 _handlers 的默认 collect() 实现，子类在 new() 中注册 (key, handler) 元组，
父类统一遍历执行并做 safe.block 容错。NvidiaGPU 采用此模式，因为 pynvml 是进程内 Python API，
每个指标可独立封装为 lambda 闭包。

其余加速器（AMD、昇腾、寒武纪等）覆写 collect() 自行实现采集逻辑，未使用 _handlers 机制。
原因是这些厂商的 CLI 工具（如 npu-smi、cnmon、hy-smi 等）单次调用返回多维度指标，
拆分为独立 lambda 会增加不必要的进程调用开销，且不利于批量解析。
覆写的 collect() 方法整体包裹在 safe.block 中，异常时返回空列表，与父类默认实现的安全语义一致。
"""

from .accelerator import ACCELERATOR_REGISTRY
from .apple import Apple
from .cpu import CPU
from .memory import Memory

__all__ = ["CPU", "Memory", "Apple", "ACCELERATOR_REGISTRY"]
