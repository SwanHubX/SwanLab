"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: gRPC 边界。

Channel 是唯一 import grpc 的地方：延迟导入、channel-ready 等待、
单次调用附加 deadline，并把 grpc 状态码映射为结构化错误。
消息上限与默认超时取自 config.py，与 Go 侧保持一致。
不持有 Popen，不做任何文件系统操作。

规划导出：Channel
"""
