#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-13 19:46:18
@File: swanlab/server/module/resp.py
@IDE: vscode
@Description:
    定义响应结构体
    说是结构体，实际上是各种函数，用于返回各种响应
    结构体名称结构为：[错误描述]_[HTTP状态码]
"""
from fastapi.responses import JSONResponse as _JSONResponse

# ---------------------------------- 错误码 ----------------------------------


_SUCCESS_200 = 0
"""一切正常，成功，期望的HTTP状态码为200"""

_PARAMS_ERROR_422 = 3422
"""参数错误，通常对应着前端传输的参数无法通过校验，这在中间件中处理，期望的HTTP状态码为422"""

_NOT_FOUND_404 = 3404
"""资源不存在，通常对应着路径不存在，期望的HTTP状态码为404"""

_CONFLICT_409 = 3409
"""可能造成冲突，通常意味着试图对资源进行不合理的操作，期望的HTTP状态码为409"""

_DATA_ERROR_500 = 3500
"""服务端存储的数据格式错误，这通常意味着指定资源无法解析为期望格式，期望的HTTP状态码为500"""

_UNEXCEPTED_ERROR_500 = 3555
"""未知错误，通常对应着未知的异常，期望的HTTP状态码为500"""

# ---------------------------------- 定义响应结构体 ----------------------------------


def _ResponseBody(code: int, message: str = None, data: dict = None):
    """构造响应，返回一个字典，包含响应码，响应信息和响应数据

    Parameters
    ----------
    code : int
        响应码，0表示成功，非0表示失败
    message : str, optional
        错误信息，如果code为0，错误信息强制为success，如果code不为0，错误信息必须提供
    data : dict, optional
        响应数据，如果传入，必须为字典类型
    """
    # 如果code为0，错误信息强制为success
    message = "success" if code == 0 else message
    # 如果code不为0，错误信息必须提供
    assert code == 0 or message is not None and len(message) > 0
    # 如果传入了响应数据，必须为字典类型
    assert data is None or isinstance(data, dict)
    # 构造响应
    if data is None:
        return {
            "code": code,
            "message": message,
        }
    else:
        return {
            "code": code,
            "message": message,
            "data": data,
        }


def SUCCESS_200(data: dict):
    """成功响应"""
    return _JSONResponse(
        status_code=200,
        content=_ResponseBody(_SUCCESS_200, data=data),
    )


def PARAMS_ERROR_422(message: str = "params error"):
    """请求参数错误"""
    return _JSONResponse(
        status_code=422,
        content=_ResponseBody(_PARAMS_ERROR_422, message=message),
    )


def NOT_FOUND_404(message: str = "NotFound"):
    """资源不存在"""
    return _JSONResponse(
        status_code=404,
        content=_ResponseBody(_NOT_FOUND_404, message=message),
    )


def CONFLICT_409(message: str = "Conflict"):
    """请求会造成资源冲突"""
    return _JSONResponse(
        status_code=409,
        content=_ResponseBody(_CONFLICT_409, message=message),
    )


def DATA_ERROR_500(message: str = "data error"):
    """服务端存储的数据格式错误"""
    return _JSONResponse(
        status_code=500,
        content=_ResponseBody(_DATA_ERROR_500, message=message),
    )


def UNEXPECTED_ERROR_500(message: str = "unexpected error"):
    """未知错误"""
    return _JSONResponse(
        status_code=500,
        content=_ResponseBody(_UNEXCEPTED_ERROR_500, message=message),
    )
