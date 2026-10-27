"""
@author: Nexisato
@file: launcher.py
@time: 2026/10/27 15:40
@description: CoreLauncher：core 服务进程的本地所有权。

唯一同时持有 process、runtime 与 channel 的对象，协议无关：
start() 负责 runtime → spawn → 等待端点 → channel-ready；
close() / abort() 负责分级回收并删除 runtime，Teardown 由 CoreGo 在调用前完成。
不认识 CoreService，也不发送任何业务 RPC；进程机制本身在 process/ 下。
超时与文件名取自 config.py，不在本模块散落硬编码值。
"""
