#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:36
@File: base.py
@IDE: pycharm
@Description:
    SwanLab OpenAPI API基类
"""
import json
from datetime import datetime, timezone
from typing import Any, Callable, List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from swanlab.api import LoginInfo
from swanlab.api.auth.login import login_by_key
from swanlab.api.http import HTTP
from swanlab.api.openapi.types import ApiResponse
from swanlab.log.log import SwanLog
from swanlab.package import get_package_version

_logger: Optional[SwanLog] = None

def get_logger(log_level: str = "info") -> SwanLog:
    global _logger
    if _logger is None:
        _logger = SwanLog("swanlab.openapi", log_level)
    else:
        _logger.level = log_level
    return _logger

def handle_response(resp: requests.Response) -> ApiResponse:
    try:
        data = resp.json()
    except (json.decoder.JSONDecodeError, requests.JSONDecodeError):
        return ApiResponse[str](
            code=resp.status_code,
            errmsg="sdk decode json error",
            data=resp.text
        )

    if not isinstance(data, dict):
        return ApiResponse[Any](
            code=resp.status_code,
            errmsg="sdk decode dict error",
            data=data
        )

    code = resp.status_code
    if 200 <= code < 300:
        message = ""
    else:
        message = f"api error: {resp.reason}. Trace id: {resp.headers.get('traceid')}"
    return ApiResponse(
        code=code,
        errmsg=message,
        data=data
    )


class ApiHTTP:
    REFRESH_TIME = HTTP.REFRESH_TIME

    def __init__(self, login_info: LoginInfo):
        self.__logger = get_logger()
        self.__login_info: LoginInfo = login_info
        self.__session: requests.Session = self.__init_session()

    @property
    def username(self) -> str:
        """
        当前登录的用户名
        """
        return self.__login_info.username or ""

    @property
    def base_url(self):
        return self.__login_info.api_host

    @property
    def sid_expired_at(self):
        """
        获取sid的过期时间
        """
        return datetime.strptime(self.__login_info.expired_at, "%Y-%m-%dT%H:%M:%S.%fZ")

    def __init_session(self) -> requests.Session:
        session = requests.Session()
        session.mount(
            prefix="https://",
            adapter=HTTPAdapter(
                max_retries=Retry(
                    total=3,
                    backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504],
                    allowed_methods=(["GET", "POST", "PUT", "DELETE", "PATCH"])
                )
            )
        )
        session.headers["swanlab-sdk"] = get_package_version()
        session.cookies.update({"sid": self.__login_info.sid})
        return session

    def __before_request(self):
        if (self.sid_expired_at - datetime.now(timezone.utc).replace(tzinfo=None)).total_seconds() < self.REFRESH_TIME:
            self.__logger.debug("Refreshing sid...")
            self.__login_info = login_by_key(self.__login_info.api_key, save=False)
            self.__session.headers["cookie"] = f"sid={self.__login_info.sid}"

    def get(self, url: str, params: dict = None) -> ApiResponse:
        self.__before_request()
        resp = self.__session.get(self.base_url + url, params=params)
        return handle_response(resp)

    def post(self, url: str, data: dict = None) -> ApiResponse:
        self.__before_request()
        resp = self.__session.post(self.base_url + url, json=data)
        return handle_response(resp)


class ApiBase:
    def __init__(self, http: ApiHTTP):
        self.http: ApiHTTP = http

def fetch_paginated_api(
    api_func: Callable[..., ApiResponse],  # 分页 API 请求函数
    page_field: str = "page",
    size_field: str = "size",
    page_size: int = 10,
    *args, **kwargs
) -> ApiResponse[List]:
    """
    通用分页全量拉取函数

    Args:
        api_func (Callable): 分页 API 请求函数，应返回 ApiResponse[Pagination]
        page_field (str): 页码参数名，默认为 "page"
        size_field (str): 每页大小参数名，默认为 "size"
        page_size (int): 每页条数，默认为 10
        *args: 传递给 api_func 的位置参数
        **kwargs: 传递给 api_func 的关键字参数

    Returns:
        ApiResponse[list]: 返回所有分页数据组成的 ApiResponse
    """
    page = 1
    objs = []
    while True:
        kwargs.update({page_field: page, size_field: page_size})
        resp: ApiResponse = api_func(*args, **kwargs)
        if resp.errmsg:
            break

        objs.extend(resp.data.list)

        if len(objs) >= resp.data.total:
            break
        page += 1

    return ApiResponse[List](code=resp.code, errmsg=resp.errmsg, data=objs)
