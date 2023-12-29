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
from io import TextIOWrapper
import math
import ujson
import os
from ..utils import create_time, get_a_lock
from urllib.parse import quote
from typing import Union
from .modules import BaseType


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
        self.__target_dict["update_time"] = create_time()
        # 此处不要添加断点，断点会导致一系列文件操作出现问题
        f.truncate()
        f.seek(0)
        ujson.dump(self.__target_dict, f, indent=4, ensure_ascii=False)

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
        with open(self.__dict_path, "w", encoding="utf-8") as f:
            self.save(f, data)


class ExperimentPoxy(object):
    """实验代理类，这类有个特点是不可以重复加载，但是会派生出其他的表单
    比如当前实验下新纪录了一个tag，则会派生出一个新的表单，用于保存这个tag的数据
    本类的主要作用是派生出其他的表单
    """

    # 每__slice_size个tag数据保存为一个文件
    __slice_size = 1000

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

    def new_tag_data(self, index) -> dict:
        """创建一个新的data数据，实际上是一个字典，包含一些默认信息"""
        return {
            "index": str(index),
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

    def save_tag(self, tag: str, data: Union[float, int, BaseType], experiment_id: int, index: int, sum: int, **kwargs):
        """保存一个tag的数据

        Parameters
        ----------
        tag : str
            tag名称
        data : _type_
            tag数据
        experiment_id : int
            实验id
        index : int
            tag索引
        sum : int
            当前tag总数
        """
        if isinstance(data, BaseType):
            data.tag = quote(tag, safe="")
            data = data.get_data()
        # 创建一个新的tag数据
        new_tag_data = self.new_tag_data(index)
        # new_tag_data["experiment_id"] = experiment_id
        new_tag_data["data"] = data
        # 对于kwargs中的数据，剔除其中value为None的数据
        for key, value in kwargs.items():
            if value is not None:
                new_tag_data[key] = value
        # 存储路径
        # tag需要转译，/会导致路径错误，需转译为%2F
        save_folder = os.path.join(self.path, quote(tag, safe=""))
        if not os.path.exists(save_folder):
            os.mkdir(save_folder)

        # 优化文件分片，每__slice_size个tag数据保存为一个文件，通过sum来判断
        need_slice = (sum - 1) % self.__slice_size == 0 or sum == 1
        mu = math.ceil(sum / self.__slice_size)
        # 存储路径
        file_path = os.path.join(save_folder, str(mu * self.__slice_size) + ".json")
        # 如果需要新增分片存储
        if need_slice:  # 达到分片条件，需要在新文件中添加新的tag数据
            file = get_a_lock(file_path, mode="w+")  # 加锁
            data = self.new_tag()
            data["data"].append(new_tag_data)
            ujson.dump(data, file, ensure_ascii=False)
        else:
            # 如果不需要新增分片存储
            file = get_a_lock(file_path, mode="r+")  # 加锁
            data = ujson.load(file)
            file.truncate()
            file.seek(0)
            # 向列表中添加新tag数据
            data["data"].append(new_tag_data)
            data["update_time"] = create_time()
            ujson.dump(data, file, ensure_ascii=False)
        file.close()
