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
        config = CosConfig(
            Region=data["region"],
            SecretId=credentials['tmpSecretId'],
            SecretKey=credentials['tmpSecretKey'],
            Token=credentials['sessionToken'],
            Scheme='https'
        )
        self.__client = CosS3Client(config)

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
                EnableMD5=False,
                # 一年
                CacheControl="max-age=31536000",
            )
        except Exception as e:
            swanlog.error("Upload error: {}".format(e))

    def upload_files(self, buffers: List[MediaBuffer]) -> Dict[str, Union[bool, List]]:
        """
        批量上传文件，keys和local_paths的长度应该相等
        :param buffers: 本地文件的二进制对象集合
        """
        pool = SimpleThreadPool()
        for buffer in buffers:
            self.upload(buffer)
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
