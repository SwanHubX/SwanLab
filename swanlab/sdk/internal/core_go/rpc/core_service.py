"""
@author: Nexisato
@file: core_service.py
@time: 2026/10/27 15:40
@description: CoreService 的 RPC 封装。

构造 protobuf 请求并复用 Channel 的 deadline 与错误映射。
所有 run 级方法显式携带 session_id。

规划导出：CoreServiceRpc
"""
