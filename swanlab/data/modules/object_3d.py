# -*- coding: utf-8 -*-
"""
Author: nexisato
Date: 2024-02-21 17:52:22
FilePath: /SwanLab/swanlab/data/modules/object_3d.py
Description:
 3D Point Cloud data parsing
"""

from .base import BaseType
import numpy as np
from typing import Union, ClassVar, Set, List
from ..utils.file import get_file_hash_numpy_array, get_file_hash_path
import os
import json
import codecs


class Object3D(BaseType):
    """Object 3D class constructor

    Parameters
    ----------
    data_or_path: numpy.array
        Path to an object3d format file or numpy array of object3d.

        numpy.array: 3D point cloud data, shape (N, 3), (N, 4) or (N, 6).
        (N, 3) :  N  * (x, y, z) coordinates
        (N, 4) :  N  * (x, y, z, c) coordinates, where c in range [1, 14]
        (N, 6) :  N  * (x, y, z, r, g, b) coordinates
    """

    SUPPORTED_TYPES: ClassVar[Set[str]] = {
        "obj",
        "gltf",
        "glb",
        "babylon",
        "stl",
        "pts.json",
    }

    def __init__(
        self,
        data_or_path: Union[np.ndarray, List["Object3D"]],
        caption: str = None,
    ):
        super().__init__(data_or_path)
        self.object3d_data = None
        self.caption = self.__convert_caption(caption)

    def get_data(self):
        # 如果传入的是Object3D类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        self.__preprocess(self.value)
        # 目前只支持numpy.array格式的3D点云数据
        hash_name = get_file_hash_numpy_array(self.object3d_data)[:16]
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = (
            f"{self.caption}-step{self.step}-{hash_name}.pts.json"
            if self.caption is not None
            else f"object3d-step{self.step}-{hash_name}.pts.json"
        )
        # 如果不存在目录则创建
        if os.path.exists(save_dir) is False:
            os.makedirs(save_dir)
        save_path = os.path.join(save_dir, save_name)
        # 保存图像到指定目录
        self.__save(save_path)
        return save_name

    def __preprocess(self, data_or_path):
        """根据输入不同的输入类型进行不同处理"""
        if isinstance(data_or_path, str):
            # TODO
            raise NotImplementedError("The input type of string is not supported yet.")

        elif isinstance(data_or_path, np.ndarray):
            # 如果输入为numpy.array，那么必须是满足我们所定义的shape的3D点云数据
            if len(data_or_path.shape) != 2 or data_or_path.shape[1] not in {3, 4, 6}:
                raise TypeError(
                    """
                    The shape of the numpy array must be one of either:
                        (N, 3) :  N  * (x, y, z) coordinates
                        (N, 4) :  N  * (x, y, z, c) coordinates, where c in range [1, 14]
                        (N, 6) :  N  * (x, y, z, r, g, b) coordinates
                    """
                )
            self.object3d_data = data_or_path
        else:
            raise TypeError("Invalid data type, only support string or numpy.array.")

    def __convert_caption(self, caption):
        """将caption转换为字符串"""
        # 如果类型是字符串，则不做转换
        if isinstance(caption, str):
            caption = caption
        # 如果类型是数字，则转换为字符串
        elif isinstance(caption, (int, float)):
            caption = str(caption)
        # 如果类型是None，则转换为默认字符串
        elif caption is None:
            caption = None
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption

    def __save(self, save_path):
        """
        保存 3D点云资源文件 .pts.json 到指定路径
        """
        try:
            list_data = self.object3d_data.tolist()
            with codecs.open(save_path, "w", encoding="utf-8") as fp:
                json.dump(
                    list_data,
                    fp,
                    separators=(",", ":"),
                    sort_keys=True,
                    indent=4,
                )
        except Exception as e:
            raise ValueError(f"Could not save the 3D point cloud data to the path: {save_path}") from e

    def expect_types(self, *args, **kwargs) -> list:
        """返回支持的文件类型"""
        return ["str", "numpy.array", "io"]

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Object3D"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.object3d
