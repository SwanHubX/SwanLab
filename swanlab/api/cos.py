#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/8 19:47
@File: cos.py
@IDE: pycharm
@Description:
    tencent cos
"""
# noinspection PyPackageRequirements
from qcloud_cos import CosConfig
# noinspection PyPackageRequirements
from qcloud_cos import CosS3Client
# noinspection PyPackageRequirements
from qcloud_cos.cos_threadpool import SimpleThreadPool
from datetime import datetime, timedelta
from typing import List, Dict, Union


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
        config = CosConfig(
            Region=data["region"],
            SecretId=credentials['tmpSecretId'],
            SecretKey=credentials['tmpSecretKey'],
            Token=credentials['sessionToken'],
            Scheme='https'
        )
        self.__client = CosS3Client(config)

    def upload(self, key: str, local_path):
        """
        上传文件，需要注意的是file_path应该为unix风格而不是windows风格
        开头不能有/
        :param key: 上传到cos的文件名称
        :param local_path: 本地文件路径，一般用绝对路径
        """
        key = self.__prefix + '/' + key
        self.__client.upload_file(
            Bucket=self.__bucket,
            Key=key,
            LocalFilePath=local_path,
            EnableMD5=False,
            progress_callback=None
        )

    def upload_files(self, keys: List[str], local_paths: List[str]) -> Dict[str, Union[bool, List]]:
        """
        批量上传文件，keys和local_paths的长度应该相等
        :param keys: 上传到cos的文件名称集合
        :param local_paths: 本地文件路径，需用绝对路径
        """
        assert len(keys) == len(local_paths), "keys and local_paths should have the same length"
        pool = SimpleThreadPool()
        for key, local_path in zip(keys, local_paths):
            pool.add_task(self.upload, key, local_path)
        pool.wait_completion()
        result = pool.get_result()
        return result

    @property
    def should_refresh(self):
        # cos传递的是北京时间，需要添加8小时
        now = datetime.utcnow() + timedelta(hours=8)
        # 过期时间减去当前时间小于刷新时间，需要注意为负数的情况
        if self.__expired_time < now:
            return True
        return (self.__expired_time - now).seconds < self.REFRESH_TIME
