"""
@author: cunyue
@file: service.py
@time: 2025/12/11 18:48
@description: 服务相关API接口
"""

from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from typing import List, Tuple

import requests

from ..client import Client
from ...log import swanlog
from ...toolkit.models.data import MediaBuffer


def _upload(*, url: str, buffer: BytesIO):
    """
    上传文件到COS
    :param url: COS上传URL
    :param buffer: 文件内容的BytesIO对象
    """
    buffer.seek(0)
    response = requests.put(url, data=buffer, headers={'Content-Type': 'application/octet-stream'})
    response.raise_for_status()


def upload_to_cos(client: Client, *, cuid: str, buffers: List[MediaBuffer]):
    """
    上传文件到COS
    :param client: 对应的客户端实例
    :param cuid: 实验cuid
    :param buffers: 媒体数据缓冲区
    """
    failed_buffers: List[Tuple[str, MediaBuffer]] = []
    # 1. 后端签名
    data, _ = client.post(
        '/resources/presigned/put',
        {"experimentId": cuid, "paths": [buffer.file_name for buffer in buffers]},
    )
    urls: List[str] = data['urls']
    # 2. 并发上传
    # executor.submit可能会失败，因为线程数有限或者线程池已经关闭
    # 来自此issue: https://github.com/SwanHubX/SwanLab/issues/889，此时需要一个个发送
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for index, buffer in enumerate(buffers):
            url = urls[index]
            try:
                # 指针回到开头
                futures.append(executor.submit(_upload, url=url, buffer=buffer))
            except RuntimeError:
                failed_buffers.append((url, buffer))
        for future in futures:
            future.result()
    # 重试失败的buffer
    if len(failed_buffers):
        swanlog.debug("Retrying failed buffers: {}".format(len(failed_buffers)))
        for url, buffer in failed_buffers:
            _upload(url=url, buffer=buffer)
