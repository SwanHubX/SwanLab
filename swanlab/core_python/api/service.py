"""
@author: cunyue
@file: service.py
@time: 2025/12/11 18:48
@description: 服务相关API接口
"""

import time
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from typing import List, Tuple

import requests
from requests.exceptions import RequestException

from ..client import Client
from ...log import swanlog
from ...toolkit.models.data import MediaBuffer


def upload_file(*, url: str, buffer: BytesIO, max_retries=3):
    """
    上传文件到COS
    :param url: COS上传URL
    :param buffer: 文件内容的BytesIO对象
    :param max_retries: 最大重试次数
    """
    # 这里也可以创建一个 Session 对象复用 TCP 连接
    with requests.Session() as session:
        for attempt in range(1, max_retries + 1):
            try:
                buffer.seek(0)
                response = session.put(
                    url,
                    data=buffer,
                    headers={'Content-Type': 'application/octet-stream'},
                    timeout=30,
                )
                response.raise_for_status()
                return
            except RequestException:
                swanlog.warning("Upload attempt {} failed for URL: {}".format(attempt, url))
                # 如果是最后一次尝试，抛出异常
                if attempt == max_retries:
                    raise
                # 简单的指数退避（等待 1s, 2s, 4s...）
                time.sleep(2 ** (attempt - 1))


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
        assert len(urls) == len(buffers), "URLs and buffers length mismatch"
        # 2.1 在线程中并发上传
        for index, buffer in enumerate(buffers):
            url = urls[index]
            try:
                future = executor.submit(upload_file, url=url, buffer=buffer)
                futures.append((future, url, buffer))
            except RuntimeError:
                failed_buffers.append((url, buffer))
        # 2.2 收集结果
        for future, url, buffer in futures:
            try:
                future.result()
            except Exception as e:
                swanlog.warning(f"Failed to upload {url}: {e}, will retry...")
                failed_buffers.append((url, buffer))
    # 3. 重试失败的buffer，重新上传
    if len(failed_buffers):
        swanlog.debug("Retrying failed buffers: {}".format(len(failed_buffers)))
        for url, buffer in failed_buffers:
            try:
                upload_file(url=url, buffer=buffer)
            except Exception as e:
                swanlog.error(f"Failed to upload {url}: {e}")
