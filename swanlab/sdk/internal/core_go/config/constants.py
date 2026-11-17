"""
@author: Nexisato
@file: constants.py
@time: 2026/10/27 15:40
@description: core_go 客户端与本地守护进程的跨模块协议常量及超时预算。

只收敛跨模块或与 Go core（portinfo / swanlab-core CLI / gRPC）对齐的契约值；
单文件内部的私有实现常量（如 Win32 API 标志位）保留在对应模块内。
"""

__all__ = [
    "BINARY_NAME",
    "CORE_BINARY_ENV",
    "RUNTIME_DIR_PREFIX",
    "RUNTIME_DIR_MODE",
    "PORT_FILE_NAME",
    "LOCK_FILE_NAME",
    "LOCK_FILE_MODE",
    "SOCKET_NAME",
    "STDERR_FILE_NAME",
    "PORT_FILE_MAX_BYTES",
    "MAX_UNIX_PATH_BYTES",
    "PORT_FILE_EOF",
    "EXIT_USAGE_ERROR",
    "EXIT_ENDPOINT_BUSY",
    "PORT_FILE_TIMEOUT",
    "ENDPOINT_POLL_INTERVAL",
    "ENDPOINT_PROBE_TIMEOUT",
    "START_LOCK_TIMEOUT",
    "START_LOCK_POLL_INTERVAL",
    "CHANNEL_READY_TIMEOUT",
    "RPC_TIMEOUT",
    "RUN_START_TIMEOUT",
    "UPSERT_TIMEOUT",
    "RUN_FINISH_TIMEOUT",
    "STATS_TIMEOUT",
    "CONFIRM_TIMEOUT",
    "DETACH_IDLE_TIMEOUT",
    "TEARDOWN_TIMEOUT",
    "WAIT_TIMEOUT",
    "TERMINATE_TIMEOUT",
    "KILL_TIMEOUT",
    "GRPC_MAX_MESSAGE_BYTES",
    "STDERR_TAIL_BYTES",
]

# ---------------------------------- 二进制与环境变量 ----------------------------------

# 发布产物名，实际路径由 runtime.locate_core_binary 定位。
BINARY_NAME = "swanlab-core"
# 覆盖包内二进制的环境变量；空值视为未设置。
CORE_BINARY_ENV = "SWANLAB_CORE_BINARY"

# ---------------------------------- Go portinfo / CLI 契约 ----------------------------------

# 私有 runtime 目录前缀与权限（0700）。
RUNTIME_DIR_PREFIX = "swanlab-core-"
RUNTIME_DIR_MODE = 0o700
# port-file / 启动锁 / UDS / 启动期 stderr 文件名，与 Go portinfo 对齐。
PORT_FILE_NAME = "core.port"
LOCK_FILE_NAME = "core.lock"
LOCK_FILE_MODE = 0o600
SOCKET_NAME = "core.sock"
STDERR_FILE_NAME = "core.stderr.log"

# port-file 文本协议边界与长度上限，对齐 Go portinfo。
PORT_FILE_MAX_BYTES = 4096
MAX_UNIX_PATH_BYTES = 103
PORT_FILE_EOF = b"EOF"

# swanlab-core CLI 退出码契约，对齐 core/cmd/swanlab-core/main.go。
EXIT_USAGE_ERROR = 2
EXIT_ENDPOINT_BUSY = 3

# ---------------------------------- 握手、超时与回收预算 ----------------------------------

# 等待 port-file 出现并解析成功的总预算。
PORT_FILE_TIMEOUT = 30.0
# 轮询 port-file 与子进程状态的间隔。
ENDPOINT_POLL_INTERVAL = 0.05
# 探测已有端点是否存活的短超时，对齐 Go endpointProbeTimeout。
ENDPOINT_PROBE_TIMEOUT = 0.2
# 等待启动锁的总预算与轮询间隔，对齐 Go portinfo 常量。
START_LOCK_TIMEOUT = 35.0
START_LOCK_POLL_INTERVAL = 0.02
# gRPC channel 进入 READY 的预算。
CHANNEL_READY_TIMEOUT = 10.0
# 通用单次 gRPC 调用预算，仅用于没有专属语义的超时。
RPC_TIMEOUT = 10.0
# DeliverRunStart：要等云端认证、项目与实验创建，预算放宽。
RUN_START_TIMEOUT = 30.0
# Upsert：数据面单批转发，失败交给 consumer 重试。
UPSERT_TIMEOUT = 5.0
# DeliverRunFinish：只提交终态，不等待上传排空。
RUN_FINISH_TIMEOUT = 15.0
# GetOperationStats：finish 阶段高频轮询，超时必须短。
STATS_TIMEOUT = 2.0
# ConfirmRunFinish：core 内部要排空上传队列，预算高于启动。
CONFIRM_TIMEOUT = 45.0
# Attach Mode 拉起 detached core 的默认空闲回收超时（秒）。
DETACH_IDLE_TIMEOUT = 300.0
# TeardownService 的预算，只发请求不等进程退出。
TEARDOWN_TIMEOUT = 5.0
# 回收子进程：先等优雅退出，超时依次升级为 terminate / kill。
WAIT_TIMEOUT = 10.0
TERMINATE_TIMEOUT = 3.0
KILL_TIMEOUT = 2.0

# ---------------------------------- gRPC 与诊断 ----------------------------------

# 单条消息上限，必须与 Go 侧 maxRecvMsgSize 保持一致。
GRPC_MAX_MESSAGE_BYTES = 4 << 20
# stderr 尾部最多读取的字节数。
STDERR_TAIL_BYTES = 4096
