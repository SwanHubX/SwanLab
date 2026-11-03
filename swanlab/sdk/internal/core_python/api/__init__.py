"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 18:19
@description: SwanLab 运行时 API 封装

绝大多数API使用 Client 对象，少部分API使用requests库直接调用
我们以rpc风格封装API，方便调用
"""

# 启动握手阶段（创建项目/开启实验）的重试次数。
# 保留 1 次重试覆盖瞬时抖动，同时避免 Cloud 不可用时 init 长时间阻塞。
STARTUP_RETRIES = 1
