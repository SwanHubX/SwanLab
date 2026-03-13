"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 14:00
@description: SwanLab SDK 内部诊断日志模块

与 console 模块互补：
- console：面向用户的终端美化输出（rich）
- log：面向开发者/运维的文件持久化诊断日志（logging）

同时 console 模块在 SWANLAB_DEBUG=true 时会将终端输出同步到诊断日志文件。

设计要点：
1. 基于 Python 标准库 logging，零额外依赖
2. 使用独立命名空间 "swanlab.internal"，不污染用户的 logging 配置
3. 延迟绑定（Lazy Attach）：log_dir 就绪前日志缓冲在内存中，就绪后一次性 flush 到文件
4. RotatingFileHandler 自动轮转，防止日志文件无限膨胀
"""

import logging
import os
from logging.handlers import MemoryHandler, RotatingFileHandler
from pathlib import Path
from typing import Optional

__all__ = ["debug", "info", "warning", "error", "critical", "bindfile", "reset"]

# ---------------------------------------------------------------------------
# Logger 初始化
# ---------------------------------------------------------------------------
_LOGGER_NAME = "swanlab.internal"
_logger = logging.getLogger(_LOGGER_NAME)
# 阻止日志向上传播到 root logger，避免干扰用户的 logging 配置
_logger.propagate = False
# 内部诊断日志默认捕获 DEBUG 及以上级别的所有信息
_logger.setLevel(logging.DEBUG)

# 日志格式：由调用方（console 模块）负责按 loguru 风格预格式化后传入，此处直接写原始消息
_FORMATTER = logging.Formatter(fmt="%(message)s")

# ---------------------------------------------------------------------------
# 内存缓冲 Handler（阶段 1：log_dir 就绪前）
# capacity 设为 1024 条，足够覆盖 login 阶段的少量请求
# flushLevel 设为 CRITICAL+1（即永远不会自动 flush），完全由 bindfile() 手动触发
# ---------------------------------------------------------------------------
_memory_handler: Optional[MemoryHandler] = None
_file_handler: Optional[RotatingFileHandler] = None
_bound = False

# ---------------------------------------------------------------------------
# 轮转配置常量
# ---------------------------------------------------------------------------
_MAX_BYTES = 10 * 1024 * 1024  # 单文件上限 10 MB
_BACKUP_COUNT = 3  # 保留 3 个备份（共 ~40 MB）
_LOG_FILENAME = "debug.log"


class SecureRotatingFileHandler(RotatingFileHandler):
    """
    安全的日志轮转处理器。
    确保每次创建/打开日志文件时，文件权限都被强制设置为 0600 (仅属主可读写)。
    主要用于防护在 Linux/macOS 共享算力集群上的敏感凭证泄露。
    """

    def _open(self):
        stream = super()._open()

        # 仅在类 Unix 系统（Linux/macOS）上应用 POSIX 权限
        if os.name == "posix":
            os.chmod(self.baseFilename, 0o600)

        # 注意：Windows 环境下多为个人开发机，且 os.chmod 无法操作 ACL，故跳过
        return stream


def reset() -> None:
    """
    重置日志模块状态，恢复到初始的内存缓冲模式。
    释放现有的文件句柄，清除所有已绑定的 Handler。
    """
    global _memory_handler, _file_handler, _bound

    # 1. 安全关闭并移除所有现有的 Handlers（防止文件句柄泄露）
    for handler in _logger.handlers[:]:
        handler.close()
        _logger.removeHandler(handler)

    # 2. 重新初始化内存缓冲 Handler
    _memory_handler = MemoryHandler(
        capacity=1024,
        flushLevel=logging.CRITICAL + 1,
    )
    _memory_handler.setFormatter(_FORMATTER)
    _logger.addHandler(_memory_handler)

    # 3. 重置状态标志
    _file_handler = None
    _bound = False


# 模块首次导入时，执行一次初始化
reset()


def bindfile(log_dir: Path) -> None:
    """
    将诊断日志绑定到文件。应在 log_dir 目录创建后调用。

    执行流程：
    1. 创建 SecureRotatingFileHandler 指向 log_dir/debug.log
    2. 将内存缓冲中的日志 flush 到文件
    3. 移除内存缓冲 Handler，后续日志直接写文件

    此方法可安全地重复调用（幂等），第二次调用会被忽略。

    :param log_dir: 日志目录路径，调用方需确保目录已创建
    :raises FileNotFoundError: 如果 log_dir 不存在
    """
    global _memory_handler, _file_handler, _bound

    if _bound:
        return

    log_dir = Path(log_dir)
    if not log_dir.exists():
        raise FileNotFoundError(f"Log directory does not exist: {log_dir}")

    # 1. 创建文件 Handler
    log_path = log_dir / _LOG_FILENAME
    _file_handler = SecureRotatingFileHandler(
        filename=str(log_path),
        maxBytes=_MAX_BYTES,
        backupCount=_BACKUP_COUNT,
        encoding="utf-8",
    )
    _file_handler.setFormatter(_FORMATTER)
    _file_handler.setLevel(logging.DEBUG)

    # 2. 将内存缓冲的日志 flush 到文件
    if _memory_handler is not None:
        _memory_handler.setTarget(_file_handler)
        _memory_handler.flush()
        # 移除内存 Handler，释放缓冲区
        _logger.removeHandler(_memory_handler)
        _memory_handler.close()
        _memory_handler = None

    # 3. 挂载文件 Handler
    _logger.addHandler(_file_handler)
    _bound = True

    _logger.debug("Diagnostic log bound to file: %s", log_path)


# ---------------------------------------------------------------------------
# 模块级日志函数
# 与 console 模块保持相同的 API 风格，方便调用方统一使用
# ---------------------------------------------------------------------------


def debug(msg: str, *args) -> None:
    """记录调试级别的诊断日志"""
    _logger.debug(msg, *args)


def info(msg: str, *args) -> None:
    """记录信息级别的诊断日志"""
    _logger.info(msg, *args)


def warning(msg: str, *args) -> None:
    """记录警告级别的诊断日志"""
    _logger.warning(msg, *args)


def error(msg: str, *args) -> None:
    """记录错误级别的诊断日志"""
    _logger.error(msg, *args)


def critical(msg: str, *args) -> None:
    """记录致命级别的诊断日志"""
    _logger.critical(msg, *args)
