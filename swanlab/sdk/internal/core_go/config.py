"""
@author: Nexisato
@file: config.py
@time: 2026/10/27 15:40
@description: core 客户端的协议级常量。

只放跨模块、或定义本地协议契约的值：
- 启动与发现：二进制名、runtime 目录前缀、port-file 名
- 握手与清理：等端点、channel-ready、单次 RPC、teardown/wait/terminate/kill
- gRPC：单条消息上限，必须与 Go 侧 maxRecvMsgSize 对齐
- 诊断：临时 stderr 文件名与有界尾部长度

阶段内部值（例如单次轮询间隔的微调）留在使用它的模块里，不在这里堆。
超时是 P0 的固定预算，不读取环境变量或用户设置。
"""
