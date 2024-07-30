#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/17 19:30
@File: task.py
@IDE: pycharm
@Description:
    打包、上传、开启任务
"""
import click
from .utils import login_init_sid, UseTaskHttp
# noinspection PyPackageRequirements
from qcloud_cos import CosConfig, CosS3Client
from swanlab.error import ApiError
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
    help="The entry file of the task, default by main.py",
)
@click.option(
    "-y",
    is_flag=True,
    help="Skip the confirmation prompt and proceed with the task launch",
)
@click.option(
    "--python",
    default="python3.10",
    nargs=1,
    type=click.Choice(["python3.8", "python3.9", "python3.10"]),
    help="The python version of the task, default by python3.10",
)
@click.option(
    "--combo",
    "-c",
    default=None,
    nargs=1,
    type=str,
    help="The plan of the task. Swanlab will use the default plan if not specified. "
         "You can check the plans in the official documentation.",
)
@click.option(
    "--name",
    "-n",
    default="Task_{}".format(datetime.now().strftime("%b%d_%H-%M-%S")),
    nargs=1,
    type=str,
    help="The name of the task, default by Task_{current_time}",
)
def launch(path: str, entry: str, python: str, name: str, combo: str, y: bool):
    """
    Launch a task!
    """
    if not entry.startswith(path):
        raise ValueError(f"Error: Entry file '{entry}' must be in directory '{path}'")
    entry = os.path.relpath(entry, path)
    # 获取访问凭证，生成http会话对象
    login_info = login_init_sid()
    print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")
    # 确认
    if not y:
        swanlog.info("Please confirm the following information:")
        swanlog.info(f"The target folder {FONT.yellow(path)} will be packaged and uploaded")
        swanlog.info(f"You have specified {FONT.yellow(entry)} as the task entry point. ")
        combo and swanlog.info(f"The task will use the combo {FONT.yellow(combo)}")
        ok = click.confirm(FONT.swanlab("Do you wish to proceed?"), abort=False)
        if not ok:
            return
    # 压缩文件夹
    memory_file = zip_folder(path)
    # 上传文件
    src = upload_memory_file(memory_file)
    # 发布任务
    ctm = CreateTaskModel(login_info.username, src, login_info.api_key, python, name, entry, combo)
    ctm.create()
    swanlog.info(f"Task launched successfully. You can use {FONT.yellow('swanlab task list')} to view the task.")


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
                arc_name = os.path.relpath(file_path.__str__(), start=dirpath)
                # 将文件添加到压缩文件中
                z.write(file_path.__str__(), arc_name)
    memory_file.seek(0)
    return memory_file


class CosClientForTask:
    def __init__(self, sts):
        region = sts["region"]
        self.bucket = sts["bucket"]
        token = sts["credentials"]["sessionToken"]
        secret_id = sts["credentials"]["tmpSecretId"]
        secret_key = sts["credentials"]["tmpSecretKey"]
        config = CosConfig(Region=region, SecretId=secret_id, SecretKey=secret_key, Token=token, Scheme="https")
        self.client = CosS3Client(config)
        self.key = sts["prefix"] + "/tasks/" + f"{int(time.time() * 1000)}.zip"

    def upload(self, buffer: io.BytesIO):
        return self.client.upload_file_from_buffer(
            Bucket=self.bucket,
            Key=self.key,
            Body=buffer,
            MAXThread=5,
            MaxBufferSize=5,
            PartSize=1
        )


class TaskBytesIO(io.BytesIO):

    def __init__(self, read_callback, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read_callback = read_callback

    def read(self, *args):
        self.read_callback(*args)
        return super().read(*args)


class TaskProgressBar:
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
    with UseTaskHttp() as http:
        sts = http.get("/user/codes/sts")
    cos = CosClientForTask(sts)
    val = memory_file.getvalue()
    progress = TaskProgressBar(len(val))
    buffer = TaskBytesIO(progress.update, val)
    t = threading.Thread(target=progress.start)
    t.start()
    cos.upload(buffer)
    t.join()
    return cos.key


class CreateTaskModel:
    def __init__(self, username, src, key, python, name, index, combo):
        """
        :param username: 用户username
        :param key: 用户的api_key
        :param src: 任务zip文件路径
        :param python: 任务的python版本
        :param name: 任务名称
        :param index: 任务入口文件
        :param combo: 任务的套餐类型
        """
        self.username = username
        self.src = src
        self.key = key
        self.python = python
        self.name = name
        self.index = index
        self.combo = combo

    def __dict__(self):
        if self.combo is not None:
            return {
                "username": self.username,
                "src": self.src,
                "index": self.index,
                "python": self.python,
                "conf": {"key": self.key},
                "combo": self.combo,
                "name": self.name
            }
        return {
            "username": self.username,
            "src": self.src,
            "index": self.index,
            "python": self.python,
            "conf": {"key": self.key},
            "name": self.name
        }

    def create(self):
        """
        创建任务
        """
        try:
            with UseTaskHttp() as http:
                http.post("/task", self.__dict__())
        except ApiError as e:
            if e.resp.status_code == 406:
                raise click.exceptions.UsageError("You have reached the maximum number of tasks")
            elif e.resp.status_code == 400:
                raise click.exceptions.UsageError("Incorrect combo name, please check it")
            raise e
