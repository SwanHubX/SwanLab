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

__all__ = [
    "BINARY_NAME",
    "PORT_FILE_NAME",
    "RUNTIME_DIR_PREFIX",
    "STDERR_FILE_NAME",
    "PORT_FILE_MAX_BYTES",
    "MAX_UNIX_PATH_BYTES",
    "RUNTIME_DIR_MODE",
    "PORT_FILE_TIMEOUT",
    "ENDPOINT_POLL_INTERVAL",
    "CHANNEL_READY_TIMEOUT",
    "RPC_TIMEOUT",
    "TEARDOWN_TIMEOUT",
    "WAIT_TIMEOUT",
    "TERMINATE_TIMEOUT",
    "KILL_TIMEOUT",
    "STDERR_TAIL_BYTES",
    "GRPC_MAX_MESSAGE_BYTES",
]

# ---------------------------------- 启动与发现 ----------------------------------

# 发布产物名，实际路径由调用方定位后传入。
BINARY_NAME = "swanlab-core"
# port-file 文件名，与 Go 侧约定一致。
PORT_FILE_NAME = "core.port"
# 私有 runtime 目录前缀，清理时据此判断归属。
RUNTIME_DIR_PREFIX = "swanlab-core-"
# 启动期临时 stderr 文件名，仅用于失败诊断。
STDERR_FILE_NAME = "core.stderr.log"

# port-file 内容长度上限，对齐 Go portinfo.MaxFileSize。
PORT_FILE_MAX_BYTES = 4096
# unix socket 路径字节上限，对齐 Go portinfo.maxUnixPathLen。
MAX_UNIX_PATH_BYTES = 103
# 私有 runtime 目录权限。
RUNTIME_DIR_MODE = 0o700

# ---------------------------------- 握手与清理 ----------------------------------

# 等待 port-file 出现并解析成功的总预算。
PORT_FILE_TIMEOUT = 30.0
# 轮询 port-file 与子进程状态的间隔。
ENDPOINT_POLL_INTERVAL = 0.05
# gRPC channel 进入 READY 的预算。
CHANNEL_READY_TIMEOUT = 10.0
# 单次 gRPC 调用预算。
RPC_TIMEOUT = 10.0
# TeardownService 的预算，只发请求不等进程退出。
TEARDOWN_TIMEOUT = 5.0
# 回收子进程：先等优雅退出。
WAIT_TIMEOUT = 10.0
# 回收子进程：超时后 terminate。
TERMINATE_TIMEOUT = 3.0
# 回收子进程：再超时后 kill。
KILL_TIMEOUT = 2.0

# ---------------------------------- gRPC ----------------------------------

# 单条消息上限，必须与 Go 侧 maxRecvMsgSize 保持一致。
GRPC_MAX_MESSAGE_BYTES = 4 << 20

# ---------------------------------- 诊断 ----------------------------------

# stderr 尾部最多读取的字节数。
STDERR_TAIL_BYTES = 4096
