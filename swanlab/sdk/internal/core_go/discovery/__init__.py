"""
@author: Nexisato
@file: __init__.py
@time: 2026/10/27 15:40
@description: 端点发现。

解析 port-file v1，产出 Endpoint（unix / tcp 二选一）并校验 pid 与目录归属。
纯函数，不建立 channel，也不启动进程。

规划导出：Endpoint / read_port_file / validate_ownership
"""
