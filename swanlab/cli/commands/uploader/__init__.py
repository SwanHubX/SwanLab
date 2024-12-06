#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/26 15:02
@File: __init__.py.py
@IDE: pycharm
@Description:
    上传模块
"""
import os
import pathlib
import sys
import threading
import time

import click
# noinspection PyPackageRequirements
from qcloud_cos.cos_exception import CosClientError, CosServiceError
# noinspection PyPackageRequirements
from qcloud_cos.cos_threadpool import SimpleThreadPool
from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
)
from swankit.log import FONT

from swanlab.cli.utils import login_init_sid, UseTaskHttp, CosUploader


class FolderProgress:
    def __init__(self, total_size: int):
        self.progress = Progress(
            TextColumn("{task.description}", justify="left"),
            BarColumn(),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "•",
            TimeRemainingColumn(),
        )
        self.current = 0
        self.total_size = total_size
        self.running = True

    def start(self, description: str):
        with self.progress as progress:
            for i in progress.track(range(self.total_size), description=description):
                if not self.running:
                    break
                if self.current > i:
                    continue
                time.sleep(0.5)
                while True:
                    if self.current > i or not self.running:
                        break

    def increase(self):
        self.current += 1

    def stop(self):
        self.running = False


class UploadFolderHandler:
    def __init__(self, uploader: CosUploader, progress: FolderProgress, retry=10):
        self.uploader = uploader
        self.progress = progress
        self.retry = retry
        self.error = None

    def __call__(self, path: str, key: str):
        if self.error is not None:
            return
        if self.uploader.should_refresh:
            self.uploader.refresh()
        error = None
        for i in range(self.retry):
            try:
                self.uploader.client.upload_file(
                    Bucket=self.uploader.bucket,
                    Key=key,
                    LocalFilePath=path,
                )
                return self.progress.increase()
            except CosClientError or CosServiceError as e:
                error = e
                continue
            except Exception as e:
                error = e
                break
        if error is not None:
            self.error = error
            raise error


@click.command()
@click.argument("path", type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True))
@click.option("--name", "-n", type=str, help="Name of the dataset to be uploaded.")
@click.option("--desc", "-d", type=str, help="Description of the dataset to be uploaded.")
def upload(path, name: str, desc: str):
    """
    Upload your 'dataset' to the cloud to
    accelerate task initiation speed.
    """
    path = os.path.abspath(path)
    name = name or os.path.basename(path)
    login_init_sid()

    # 创建数据集索引
    with UseTaskHttp() as http:
        data = {"name": name, "desc": desc} if desc else {"name": name}
        dataset = http.post("/task/dataset", data=data)

    # 创建上传对象
    uploader = CosUploader()
    cuid = dataset['cuid']
    prefix = f"{uploader.prefix}/datasets/{cuid}/{name}"

    # 上传文件夹
    if os.path.isdir(path):
        total_size = 0
        for root, dirs, files in os.walk(path):
            total_size += len(files)
        progress = FolderProgress(total_size)
        t = threading.Thread(target=progress.start, args=("Uploading...",))
        t.start()
        pool = SimpleThreadPool(5)
        handler = UploadFolderHandler(uploader, progress)
        # 遍历，添加任务
        for root, dirs, files in os.walk(path):
            for file in files:
                local_path = os.path.join(root, file).__str__()
                # 生成key
                key = prefix
                tmp = os.path.relpath(local_path, start=path)
                path_obj = pathlib.Path(tmp.__str__())
                folders = [parent for parent in path_obj.parents]
                folders.reverse()
                for folder in folders[1:]:
                    key += "/" + folder.name
                key += "/" + path_obj.name
                pool.add_task(handler, local_path, key)
        pool.wait_completion()
        result = pool.get_result()
        progress.stop()
        t.join()
        if not result['success_all']:
            print("Not all files upload success. you should retry, Error: {}".format(handler.error), file=sys.stderr)
            status = "FAILURE"
        else:
            status = "SUCCESS"
            print(FONT.swanlab("Upload success, dataset id: {}".format(FONT.bold(cuid))))
        return http.patch("/task/dataset/status", {"cuid": cuid, "status": status})
