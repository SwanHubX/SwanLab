#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 01:25:31
@File: swanlab\server\api\utils.py
@IDE: vscode
@Description:
    工具文件，重复使用的逻辑的封装
"""
from typing import Callable
from fastapi import Response
import ujson


def ResponseBody(code: int, message: str = None, data: dict = None):
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


async def get_response_body(response: Response, callback: Callable[[dict], dict] = None) -> Response:
    """从响应对象中获取响应体，并转换为字典，进行回调处理

    Parameters
    ----------
    response : Response
        响应对象
    callback : Callable[[dict], dict], optional
        回调函数，由于本函数本身异步，建议传入的回调函数是同步的, by default None

    Returns
    -------
    Response
        新的响应对象
    """
    response_body = b""
    async for chunk in response.body_iterator:
        response_body += chunk
    body = ujson.loads(response_body.decode("utf-8"))
    if callback is not None:
        body = callback(body)
    return Response(
        content=ujson.dumps(body).encode("utf-8"),
        status_code=response.status_code,
        headers=dict(response.headers),
        media_type=response.media_type,
    )
