"""
@author: Nexisato
@file: core.py
@time: 2026/10/27 15:40
@description: CoreGo：CoreProtocol 在 Go core 上的适配层。

只做 run 级编排：start 时 spinup + DeliverRunStart 并保存 session_id，
upsert 直接转发到 CoreServiceRpc，finish / stats / confirm 一一映射。
不做本地存储、指标上下文、批量与重试——那些是 Go core 或上游消费线程的职责。

规划导出：CoreGo
"""
