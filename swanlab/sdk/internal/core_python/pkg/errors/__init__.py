"""
@author: likang
@file: __init__.py
@time: 2026/11/3
@description: CorePython 启动握手的错误分类

该模块只服务于 Python Core 实现。跨进程 / gRPC 边界只能传输 CoreErrorKind，
不能传输异常对象。
"""

from requests.exceptions import InvalidSchema, InvalidURL, MissingSchema, RequestException

from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    CORE_ERROR_KIND_AUTH,
    CORE_ERROR_KIND_OTHER,
    CORE_ERROR_KIND_TRANSIENT,
    CoreErrorKind,
)

__all__ = ["classify_start_error"]

# URL / 协议配置错误，重试或降级都无法解决
_CONFIG_ERRORS = (InvalidURL, MissingSchema, InvalidSchema)


def classify_start_error(exc: BaseException) -> CoreErrorKind:
    """将启动握手的底层异常转换为可序列化的错误分类。"""
    if isinstance(exc, AuthenticationError):
        return CORE_ERROR_KIND_AUTH
    if isinstance(exc, ApiError):
        status = getattr(exc.response, "status_code", None)
        if status in (401, 403):
            return CORE_ERROR_KIND_AUTH
        # 没有响应通常意味着请求未完成或响应无法解析，按服务异常处理
        if status is None or status in (408, 429) or status >= 500:
            return CORE_ERROR_KIND_TRANSIENT
        return CORE_ERROR_KIND_OTHER
    # 注意：InvalidURL 等配置错误是 RequestException 的子类，必须先判断，
    # 避免坏 URL 被当成瞬时错误降级。
    if isinstance(exc, _CONFIG_ERRORS):
        return CORE_ERROR_KIND_OTHER
    # RequestException 继承自 IOError/OSError，必须先于 OSError 判断
    if isinstance(exc, RequestException):
        return CORE_ERROR_KIND_TRANSIENT
    # 本地磁盘 / 权限错误
    if isinstance(exc, OSError):
        return CORE_ERROR_KIND_OTHER
    # 业务语义错误（如 resume="must" 但实验不存在）
    if isinstance(exc, RuntimeError):
        return CORE_ERROR_KIND_OTHER
    # 200 响应缺少字段、结构异常等：decode_response 会退化为 text，随后按 dict 取值报错
    if isinstance(exc, (KeyError, TypeError)):
        return CORE_ERROR_KIND_TRANSIENT
    # 其余（ValueError / AssertionError / ValidationError 等）属于确定性的配置、
    # 后端兼容性或 SDK 内部错误，重试无意义，不参与降级
    return CORE_ERROR_KIND_OTHER
