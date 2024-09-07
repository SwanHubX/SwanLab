#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/25 14:33
@File: utils.py
@IDE: pycharm
@Description:
    一些工具函数
"""
import io
import sys
import threading
import time
from datetime import datetime, timedelta
from typing import Optional, Tuple

# noinspection PyPackageRequirements
from qcloud_cos import CosConfig, CosS3Client
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

from swanlab.api import terminal_login, create_http, LoginInfo, get_http
from swanlab.error import KeyFileError, ApiError
from swanlab.log import swanlog
from swanlab.package import get_key


class UseTaskHttp:
    """
    主要用于检测http响应是否为3xx字段，如果是则要求用户更新版本
    使用此类之前需要先调用login_init_sid()函数完成全局http对象的初始化
    """

    def __init__(self):
        self.http = get_http()

    def __enter__(self):
        return self.http

    def __exit__(self, exc_type, exc_val: Optional[ApiError], exc_tb):
        if exc_type is ApiError:
            # api已过期，需要更新swanlab版本
            if exc_val.resp.status_code // 100 == 3:
                swanlog.info("SwanLab in your environment is outdated. Upgrade: `pip install -U swanlab`")
                sys.exit(3)
        return False


def login_init_sid() -> LoginInfo:
    key = None
    try:
        key = get_key()
    except KeyFileError:
        pass
    login_info = terminal_login(key)
    create_http(login_info)
    return login_info


class UploadBytesIO(io.BytesIO):
    """
    封装BytesIO，使其可以在上传文件时显示进度条
    """

    class UploadProgressBar:
        def __init__(self, total_size: int):
            """
            :param total_size: 总大小（bytes）
            """
            self.total_size = total_size
            self.current = 0
            self.progress = Progress(
                TextColumn("{task.description}", justify="left"),
                BarColumn(),
                "[progress.percentage]{task.percentage:>3.1f}%",
                "•",
                DownloadColumn(),
                "•",
                TransferSpeedColumn(),
                "•",
                TimeRemainingColumn(),
            )

        def update(self, *args):
            self.current += args[0]

        def start(self, description: str):
            with self.progress as progress:
                for i in progress.track(range(self.total_size), description=description):
                    if self.current > i:
                        continue
                    time.sleep(0.5)
                    while True:
                        if self.current > i:
                            break

    def __init__(self, description: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.upload_progress = self.UploadProgressBar(len(self.getvalue()))
        self.t = None
        self.description = description

    def read(self, *args):
        self.upload_progress.update(*args)
        return super().read(*args)

    def __enter__(self):
        self.t = threading.Thread(target=self.upload_progress.start, args=(self.description,))
        self.t.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.t.join()
        return False


class CosUploader:
    REFRESH_TIME = 60 * 60 * 1.5  # 1.5小时

    def __init__(self):
        """初始化 cos"""
        client, sts = self.create()
        self.__expired_time = datetime.fromtimestamp(sts["expiredTime"])
        self.prefix = sts["prefix"]
        self.bucket = sts["bucket"]
        self.client = client
        self.__updating = False
        """
        标记是否正在更新sts
        """
        self.__token = sts["credentials"]["sessionToken"]  # 临时密钥使用的 token

    @property
    def token(self):
        if self.should_refresh:
            self.refresh()
        return self.__token

    @property
    def should_refresh(self):
        # cos传递的是北京时间，需要添加8小时
        now = datetime.utcnow() + timedelta(hours=8)
        # 过期时间减去当前时间小于刷新时间，需要注意为负数的情况
        if self.__expired_time < now:
            return True
        return (self.__expired_time - now).seconds < self.REFRESH_TIME

    @staticmethod
    def create() -> Tuple[CosS3Client, dict]:
        with UseTaskHttp() as http:
            sts = http.get("/user/codes/sts")
            region = sts["region"]
            token = sts["credentials"]["sessionToken"]
            secret_id = sts["credentials"]["tmpSecretId"]
            secret_key = sts["credentials"]["tmpSecretKey"]
            config = CosConfig(Region=region, SecretId=secret_id, SecretKey=secret_key, Token=token, Scheme="https")
            client = CosS3Client(config)
        return client, sts

    def refresh(self):
        """
        更新sts
        """
        # 防止多线程更新sts
        if self.__updating:
            while self.__updating:
                time.sleep(1)
            return

        self.__updating = True
        client, sts = self.create()
        self.client = client
        self.__expired_time = datetime.fromtimestamp(sts["expiredTime"])
        self.prefix = sts["prefix"]
        self.bucket = sts["bucket"]
        self.__token = sts["credentials"]["sessionToken"]
