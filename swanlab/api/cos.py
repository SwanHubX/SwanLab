#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/8 19:47
@File: cos.py
@IDE: pycharm
@Description:
    cloud object storage
"""
from concurrent.futures.thread import ThreadPoolExecutor
from datetime import datetime, timedelta
from typing import List

import boto3
from botocore.config import Config as BotocoreConfig

from swanlab.data.modules import MediaBuffer
from swanlab.log import swanlog


class CosClient:
    REFRESH_TIME = 60 * 60 * 1.5  # 1.5小时

    def __init__(self, data):
        """
        初始化cos客户端
        """
        self.__expired_time = datetime.fromtimestamp(data["expiredTime"])
        self.__prefix = data["prefix"]
        self.__bucket = data["bucket"]
        credentials = data["credentials"]

        # 往期版本适配
        end_point = data.get("endPoint", None)
        path_style = 'path' if data.get('pathStyle', False) else 'virtual'
        endpoint_url = f"https://cos.{data['region']}.myqcloud.com" if end_point is None else end_point

        self.__client = boto3.client(
            's3',
            endpoint_url=endpoint_url,
            api_version='2006-03-01',
            aws_access_key_id=credentials['tmpSecretId'],
            aws_secret_access_key=credentials['tmpSecretKey'],
            aws_session_token=credentials['sessionToken'],
            config=BotocoreConfig(signature_version="s3", s3={'addressing_style': path_style}),
        )

    def upload(self, buffer: MediaBuffer):
        """
        上传文件，需要注意的是file_path应该为unix风格而不是windows风格
        开头不能有/
        :param buffer: 本地文件的二进制数据
        """
        key = "{}/{}".format(self.__prefix, buffer.file_name)
        try:
            swanlog.debug("Uploading file: {}".format(key))
            self.__client.put_object(
                Bucket=self.__bucket,
                Key=key,
                Body=buffer.getvalue(),
                # 一年
                CacheControl="max-age=31536000",
            )
        except Exception as e:
            swanlog.error("Upload error: {}".format(e))

    def upload_files(self, buffers: List[MediaBuffer]):
        """
        批量上传文件，keys和local_paths的长度应该相等
        :param buffers: 本地文件的二进制对象集合
        """
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(self.upload, buffer) for buffer in buffers]
            for future in futures:
                future.result()

    @property
    def should_refresh(self):
        # cos传递的是北京时间，需要添加8小时
        # FIXME Use timezone-aware objects to represent datetimes in UTC; e.g. by calling .now(datetime.UTC)
        now = datetime.utcnow() + timedelta(hours=8)
        # 过期时间减去当前时间小于刷新时间，需要注意为负数的情况
        if self.__expired_time < now:
            return True
        return (self.__expired_time - now).seconds < self.REFRESH_TIME
