#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-01 20:35:49
@File: swanlab/database/table.py
@IDE: vscode
@Description:
    数据库表单类，用于操作和记录数据库表单相关的信息，实际上是一个MutableMapping，可以像字典一样操作，但是对字典数据设置的时候进行一些特殊行为
"""
from collections.abc import MutableMapping
from typing import Any
import portalocker  # 文件锁
from io import TextIOWrapper
import ujson
import os
from ..utils import create_time
from ..utils import lock_file, get_a_lock
import sys


class ProjectTablePoxy(MutableMapping):
    """项目表单代理类，这类的表单有个特点是可以重复加载，并且需要保存到文件中"""

    def __init__(self, data: dict, path: str):
        # 为data添加一些默认信息，如创建时间等
        # 判断create_time和update_time是否存在，都不存在，则添加，都存在，跳过，否则抛出异常
        if "create_time" not in data and "update_time" not in data:
            time = create_time()
            data["create_time"] = time
            data["update_time"] = time
        elif "create_time" in data and "update_time" in data:
            pass
        else:
            raise ValueError("invalid table data")
        # 保存表单信息
        self.__target_dict = data
        self.__dict_path = path

    def __getitem__(self, key):
        return self.__target_dict[key]

    def __setitem__(self, key: str, value: Any):
        """当字典数据发生变化时，需要将数据写入到文件中

        Parameters
        ----------
        key : str
            数据的名称
        value :
            数据的值，期望是一个可以被序列化的对象，但是这里不做检查
        """
        self.__target_dict[key] = value

    def __delitem__(self, key):
        del self.__target_dict[key]

    def __iter__(self):
        return iter(self.__target_dict)

    def __len__(self):
        return len(self.__target_dict)

    def __dict__(self):
        return self.__target_dict

    def save(self, f: TextIOWrapper, data: dict = None):
        """将数据保存到文件中
        是否加锁取决于传入的文件对象是否加锁，这里不做要求
        可以选择将data传入，此时会将data保存到文件中，并且更新对象信息为data信息
        """
        if data is not None:
            self.__target_dict = data
        # 此处不要添加断点，断点会导致一系列文件操作出现问题
        f.truncate()
        f.seek(0)
        ujson.dump(self.__target_dict, f, indent=4)

    def save_with_lock(self, data: dict = None):
        """将数据保存到文件中，并加锁
        可以选择将data传入，此时会将data保存到文件中，并且更新对象信息为data信息
        """
        with get_a_lock(self.__dict_path, mode="a+") as f:
            self.save(f, data)

    def save_no_lock(self, data: dict = None):
        """将数据保存到文件中，并加锁
        可以选择将data传入，此时会将data保存到文件中，并且更新对象信息为data信息
        """
        with open(self.__dict_path, "w") as f:
            self.save(f, data)


class ExperimentPoxy(object):
    """实验代理类，这类有个特点是不可以重复加载，但是会派生出其他的表单
    比如当前实验下新纪录了一个tag，则会派生出一个新的表单，用于保存这个tag的数据
    本类的主要作用是派生出其他的表单
    """

    def __init__(self, path: str):
        """初始化一个实验，保存实验信息和实验派生的文件路径

        Parameters
        ----------
        name : str
            实验名称
        path : str
        """
        self.path = path
        time = create_time()
        self.create_time = time
        self.update_time = time

    def new_tag_data(self) -> dict:
        """创建一个新的data数据，实际上是一个字典，包含一些默认信息"""
        return {
            "create_time": create_time(),
        }

    def new_tag(self) -> dict:
        """创建一个新的tag data数据集合

        Returns
        -------
        dict
            返回一个新的data数据集合
        """
        time = create_time()
        return {
            "create_time": time,
            "update_time": time,
            "data": [],
        }

    def save_tag(self, tag: str, data: Any, experiment_id: int):
        """保存一个tag的数据

        Parameters
        ----------
        tag : str
            tag名称
        data : _type_
            tag数据
        path : str
            文件路径
        """
        # 创建一个新的tag数据
        new_tag_data = self.new_tag_data()
        new_tag_data["experiment_id"] = experiment_id
        new_tag_data["data"] = data
        # 存储路径
        path = os.path.join(self.path, f"{tag}.json")
        # 拿到原本的数据
        if not os.path.exists(path):
            ujson.dump(self.new_tag(), open(path, "w"))
        # 读取数据
        with open(path, "r") as f:
            tag_data = ujson.load(f)
            tag_data["data"].append(new_tag_data)
            tag_data["update_time"] = new_tag_data["create_time"]
        # 保存数据
        with open(path, "w") as f:
            ujson.dump(tag_data, f)
