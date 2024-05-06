#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/5 22:47
@File: _utils.py
@IDE: pycharm
@Description:
    工具函数
"""
import io
import hashlib


def get_file_hash_path(file_path: str) -> str:
    """计算并返回给定文件的SHA-256哈希值。"""

    hash_sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:  # 以二进制读取模式打开文件
        while chunk := f.read(8192):  # 读取文件的小块进行处理
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()


def get_file_hash_numpy_array(array) -> str:
    """计算并返回给定NumPy数组的SHA-256哈希值。"""

    hash_sha256 = hashlib.sha256()
    # 将NumPy数组转换为字节串，然后更新哈希值
    hash_sha256.update(array.tobytes())
    return hash_sha256.hexdigest()


def get_file_hash_pil(image) -> str:
    """计算并返回给定PIL.Image对象的SHA-256哈希值。"""

    hash_sha256 = hashlib.sha256()
    # 将图像转换为字节数据
    with io.BytesIO() as buffer:
        image.save(buffer, format="PNG")  # 可以选择其他格式，如'JPEG'
        hash_sha256.update(buffer.getvalue())
    return hash_sha256.hexdigest()
