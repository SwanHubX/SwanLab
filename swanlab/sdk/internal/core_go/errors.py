"""
@author: Nexisato
@file: errors.py
@time: 2026/10/27 15:40
@description: 结构化错误与原因码。

CoreGoError 携带阶段（stage）、原因码（reason）与可选的 gRPC 状态码，
让上层能区分启动失败、通道超时、RPC 拒绝与 UNIMPLEMENTED。
内部专用，不复用 swanlab/exceptions.py 中的 HTTP/DataStore 异常。

规划导出：CoreGoError / ReasonCode
"""
