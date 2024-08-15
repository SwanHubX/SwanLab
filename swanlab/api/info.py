#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-20 14:41:25
@File: swanlab/auth/info.py
@IDE: vscode
@Description:
    定义认证数据格式
"""
import requests
from swanlab.package import save_key
from typing import Union


class LoginInfo:
    """
    登录信息类，负责解析登录接口返回的信息，并且进行保存
    无论接口请求成功还是失败，都会初始化一个LoginInfo对象
    """

    def __init__(self, resp: requests.Response, api_key: str):
        self.__api_key = api_key
        self.__resp = resp
        self.__body = resp.json() if resp.status_code == 200 else {}
        self.__username = None
        """
        如果此属性不为None，username返回此属性
        """

    @property
    def sid(self) -> Union[str, None]:
        """
        获取sid，如果请求失败则返回None
        """
        return self.__body.get("sid")

    @property
    def expired_at(self) -> Union[str, None]:
        """
        获取过期时间，如果请求失败则返回None
        """
        return self.__body.get("expiredAt")

    @property
    def username(self) -> Union[str, None]:
        """
        获取用户名，如果请求失败则返回None
        """
        if self.__username is not None:
            return self.__username
        return self.__body.get("userInfo", {}).get("username")

    @username.setter
    def username(self, value):
        """
        设置用户名
        """
        self.__username = value

    @property
    def is_fail(self):
        """
        判断登录是否失败
        """
        return self.__resp.status_code != 200

    @property
    def api_key(self):
        """
        获取api_key
        """
        if self.is_fail:
            return None
        return self.__api_key

    def __str__(self) -> str:
        """错误时会返回错误信息"""
        if self.__resp.reason == "OK":
            return "Login success"
        if self.__resp.reason == "Unauthorized" or self.__resp.reason == "Authorization Required":
            return "Error api key"
        if self.__resp.reason == "Forbidden":
            return "You need to be verified first"
        return str(self.__resp.status_code) + " " + self.__resp.reason

    def save(self):
        """
        保存登录信息
        """
        return save_key("user", self.api_key)


class ProjectInfo:
    def __init__(self, data: dict):
        self.__data = data

    @property
    def cuid(self):
        return self.__data["cuid"]

    @property
    def name(self):
        return self.__data["name"]

    @property
    def history_exp_count(self):
        return self.__data.get('_count', {'experiments': 0})["experiments"]


class ExperimentInfo:

    def __init__(self, data: dict):
        self.__data = data

    @property
    def cuid(self):
        return self.__data["cuid"]

    @property
    def name(self):
        return self.__data["name"]
