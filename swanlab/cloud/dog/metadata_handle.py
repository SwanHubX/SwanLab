#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/5 18:20
@File: metadata_handle.py
@IDE: pycharm
@Description:
    元数据处理器，看门狗嗅探元数据文件夹，向聚合器发送元数据信息
"""
from ..task_types import UploadType
from .sniffer_queue import SnifferQueue
from typing import Union, List
from watchdog.events import FileSystemEventHandler, FileSystemEvent
from swanlab.log import swanlog
import os
from queue import Queue


class MetaHandle(FileSystemEventHandler):
    ValidFiles = ['config.yaml', 'requirements.txt', 'swanlab-metadata.json']
    """
    有效的元数据文件列表，只有这些文件会被传输，如果出现其他文件出现waning
    """

    ModifiableFiles = [ValidFiles[0], ValidFiles[2]]
    """
    可修改的元数据文件列表（其他只会传输一次）
    """

    def __init__(self, queue: Queue, watched_path: str):
        """
        初始化日志嗅探处理器
        :param watched_path: 监听的路径，用作初始对照
        """
        self.watched_path = watched_path
        self.queue = SnifferQueue(queue, readable=False)
        self.on_init_upload()

    def list_all_meta_files(self) -> List[str]:
        """
        列出所有的元数据文件
        """
        files = [x for x in os.listdir(self.watched_path) if os.path.isfile(self.fmt_file_path(x)[0])]
        return [x for x in files if x in self.ValidFiles]

    def fmt_file_path(self, file_name: Union[List[str], str]) -> List[str]:
        """
        格式化文件路径
        """
        if isinstance(file_name, str):
            file_name = [file_name]
        return [os.path.join(self.watched_path, x) for x in file_name]

    def on_init_upload(self):
        """
        实例化的时候进行一次文件扫描，watched_path下所有ValidFiles生成一个一个msg发给队列
        """
        meta_files = self.list_all_meta_files()
        if len(meta_files) == 0:
            return swanlog.warning("empty meta files, it might be a bug?")
        self.queue.put((self.fmt_file_path(meta_files), UploadType.FILE))

    def on_modified(self, event: FileSystemEvent) -> None:
        """
        文件被修改时触发
        """
        if event.is_directory:
            return
        file_name = os.path.basename(event.src_path)
        if file_name not in self.ModifiableFiles:
            # 被忽略
            return swanlog.warning(f"file {file_name} is not allowed to be modified")
        self.queue.put((self.fmt_file_path(file_name), UploadType.FILE))
