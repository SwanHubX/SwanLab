"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: 子进程与 runtime 租约。

负责私有 runtime 目录、Popen 启动、等待端点、进程早死检测与
wait → terminate → kill 的分级回收。不构造 protobuf，不解释 session_id。

规划导出：RuntimeLease / acquire_runtime / spawn_core /
wait_endpoint / reap_process / cleanup_runtime
"""
