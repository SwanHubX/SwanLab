#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/17 19:30
@File: job.py
@IDE: pycharm
@Description:
    打包、上传、开启任务
"""
import click
from .utils import login_init_sid, TaskModel
from swanlab.api import get_http
# noinspection PyPackageRequirements
from qcloud_cos import CosConfig, CosS3Client
from swanlab.log import swanlog
from swankit.log import FONT
import zipfile
import threading
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)
from datetime import datetime
import time
import io
import os


@click.command()
@click.argument(
    "path",
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    default=".",
    nargs=1,
    required=True,
)
@click.option(
    "--entry",
    "-e",
    default="main.py",
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        resolve_path=True,
        readable=True,
    ),
    help="The entry file of the job, default by main.py",
)
@click.option(
    "--python",
    default="python3.10",
    nargs=1,
    type=click.Choice(["python3.8", "python3.9", "python3.10"]),
    help="The python version of the job, default by python3.10",
)
@click.option(
    "--name",
    "-n",
    default="Job_{}".format(datetime.now().strftime("%b%d_%H-%M-%S")),
    nargs=1,
    type=str,
    help="The name of the job, default by Job_{current_time}",
)
def launch(path: str, entry: str, python: str, name: str):
    if not entry.startswith(path):
        raise ValueError(f"Error: Entry file '{entry}' must be in directory '{path}'")
    entry = entry[len(path):]
    # 获取访问凭证，生成http会话对象
    login_info = login_init_sid()
    print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")
    # 上传文件
    text = f"The target folder {FONT.yellow(path)} will be packaged and uploaded, "
    text += f"and you have specified {FONT.yellow(entry)} as the task entry point. "
    swanlog.info(text)
    # click.confirm(FONT.swanlab("Do you wish to proceed?"))
    # 压缩文件夹
    memory_file = zip_folder(path)
    # 上传文件
    src = upload_memory_file(memory_file)
    # 发布任务
    http = get_http()
    http.base_url = "http://172.16.42.24"
    http.session.headers.update({'payload': str({'uid': 1, 'username': 'cunyue'})})
    get_http().post("/task", TaskModel(login_info.username, src, login_info.api_key, python, name, entry).__dict__())
    swanlog.info(f"Job launched successfully. You can use {FONT.yellow('swanlab job list')} to view the job.")


def zip_folder(dirpath: str) -> io.BytesIO:
    """
    压缩文件夹
    :param dirpath: 传入文件夹路径
    """
    memory_file = io.BytesIO()
    z = zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED)
    fs = [x for x in os.walk(dirpath)]

    progress = Progress(
        TextColumn("{task.description}", justify="left"),
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.1f}%",
        "•",
        TimeRemainingColumn(),
    )
    with progress:
        for i in progress.track(range(len(fs)), description=FONT.swanlab("Packing...  ")):
            root, dirs, files = fs[i]
            for file in files:
                # 构建文件的完整路径
                file_path = os.path.join(root, file)
                # 构建在压缩文件中的路径
                arcname = os.path.relpath(file_path.__str__(), start=dirpath)
                # 将文件添加到压缩文件中
                z.write(file_path.__str__(), arcname)
    memory_file.seek(0)
    return memory_file


class CosClientForJob:
    def __init__(self, sts):
        region = sts["region"]
        self.bucket = sts["bucket"]
        token = sts["credentials"]["sessionToken"]
        secret_id = sts["credentials"]["tmpSecretId"]
        secret_key = sts["credentials"]["tmpSecretKey"]
        config = CosConfig(Region=region, SecretId=secret_id, SecretKey=secret_key, Token=token, Scheme="https")
        self.client = CosS3Client(config)
        self.key = sts["prefix"] + "/jobs/" + f"{int(time.time() * 1000)}.zip"

    def upload(self, buffer: io.BytesIO):
        time.sleep(1)  # cos那边的密钥甚至不是立即生效，需要等一下子，否则有机率出问题
        return self.client.upload_file_from_buffer(
            Bucket=self.bucket,
            Key=self.key,
            Body=buffer,
            MAXThread=5,
            MaxBufferSize=5,
            PartSize=1
        )


class JobBytesIO(io.BytesIO):

    def __init__(self, read_callback, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read_callback = read_callback

    def read(self, *args):
        self.read_callback(*args)
        return super().read(*args)


class JobProgressBar:
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

    def start(self):
        with self.progress as progress:
            for i in progress.track(range(self.total_size), description=FONT.swanlab("Uploading...")):
                if self.current > i:
                    continue
                time.sleep(0.5)
                while True:
                    if self.current > i:
                        break


def upload_memory_file(memory_file: io.BytesIO) -> str:
    """
    上传内存文件
    :returns 上传成功后的文件路径
    """
    sts = get_http().get("/user/codes/sts")
    cos = CosClientForJob(sts)
    val = memory_file.getvalue()
    progress = JobProgressBar(len(val))
    buffer = JobBytesIO(progress.update, val)
    t = threading.Thread(target=progress.start)
    t.start()
    src = cos.upload(buffer)["Location"]
    t.join()
    return src
