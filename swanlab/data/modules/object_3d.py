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
from typing import Union, ClassVar, Set, List, Optional
from ._utils import get_file_hash_numpy_array, get_file_hash_path
import os
import json
import shutil
from io import BytesIO

# 格式化输出 json
import codecs


class Object3D(BaseType):
    """Object 3D class constructor

    Parameters
    ----------
    data_or_path: numpy.array, string, io
        Path to an object3d format file or numpy array of object3d or Bytes.IO.

        numpy.array: 3D point cloud data, shape (N, 3), (N, 4) or (N, 6).
        (N, 3) :  N  * (x, y, z) coordinates
        (N, 4) :  N  * (x, y, z, c) coordinates, where c in range [1, 14]
        (N, 6) :  N  * (x, y, z, r, g, b) coordinates
    caption: str
        caption associated with the object3d for display
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
        data_or_path: Union[str, "np.ndarray", "BytesIO", List["Object3D"]],
        caption: Optional[str] = None,
    ):
        super().__init__(data_or_path)
        self.object3d_data = None
        self.caption = self.__convert_caption(caption)
        self.extension = None

    def get_data(self):
        # 如果传入的是Object3D类列表
        if isinstance(self.value, list):
            return self.get_data_list()

        self.object3d_data = self.__preprocess(self.value)

        # 判断是否要保存(mode='disabled'时不保存)
        if not self.settings.should_save:
            return

        # 根据不同的输入类型进行不同的哈希校验
        hash_name = (
            get_file_hash_numpy_array(self.object3d_data)[:16]
            if isinstance(self.object3d_data, np.ndarray)
            else get_file_hash_path(self.object3d_data)[:16]
        )

        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = save_name = f"obj3d-step{self.step}-{hash_name}.{self.format}"
        # 如果不存在目录则创建
        if os.path.exists(save_dir) is False:
            os.makedirs(save_dir)
        save_path = os.path.join(save_dir, save_name)

        self.__save(save_path)
        return save_name

    def __preprocess(self, data_or_path):
        """根据输入不同的输入类型进行不同处理"""
        # 如果类型为 str，进行文件后缀格式检查
        if isinstance(data_or_path, str):
            extension = None
            for SUPPORTED_TYPE in Object3D.SUPPORTED_TYPES:
                if data_or_path.endswith(SUPPORTED_TYPE):
                    extension = SUPPORTED_TYPE
                    break
            if not extension:
                raise TypeError(
                    "File '"
                    + data_or_path
                    + "' is not compatible with Object3D: supported types are: "
                    + ", ".join(Object3D.SUPPORTED_TYPES)
                )
            self.extension = extension
            return data_or_path
        # 如果类型为 io.BytesIO 二进制流，直接返回
        elif isinstance(data_or_path, BytesIO):
            self.extension = "pts.json"
            return data_or_path

        # 如果类型为 numpy.array，进行numpy格式检查
        elif isinstance(data_or_path, np.ndarray):
            if len(data_or_path.shape) != 2 or data_or_path.shape[1] not in {3, 4, 6}:
                raise TypeError(
                    """
                    The shape of the numpy array must be one of either:
                        (N, 3) :  N  * (x, y, z) coordinates
                        (N, 4) :  N  * (x, y, z, c) coordinates, where c in range [1, 14]
                        (N, 6) :  N  * (x, y, z, r, g, b) coordinates
                    """
                )
            self.extension = "pts.json"
            return data_or_path
        else:
            raise TypeError("swanlab.Object3D accepts a file path or numpy like data as input")

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
        return caption.strip() if caption else None

    def __save_numpy(self, save_path):
        """保存 numpy.array 格式的 3D点云资源文件 .pts.json 到指定路径"""
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
            raise TypeError(f"Could not save the 3D point cloud data to the path: {save_path}") from e

    def __save(self, save_path):
        """
        保存 3D点云资源文件到指定路径
        """
        if isinstance(self.object3d_data, str):
            shutil.copy(self.object3d_data, save_path)
        elif isinstance(self.object3d_data, BytesIO):
            with open(save_path, "wb") as f:
                f.write(self.object3d_data.read())
        elif isinstance(self.object3d_data, np.ndarray):
            self.__save_numpy(save_path)

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        # 如果传入的是Objet3d类列表
        if isinstance(self.value, list):
            return self.get_more_list()
        else:
            return (
                {
                    "caption": self.caption,
                }
                if self.caption is not None
                else None
            )

    def expect_types(self, *args, **kwargs) -> list:
        """返回支持的文件类型"""
        return ["str", "numpy.array", "io"]

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Object3D"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.object3d
