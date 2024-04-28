#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-12 01:47:34
@File: swanlab\data\utils\file.py
@IDE: vscode
@Description:
    文件/格式检查和操作
"""
import os
import platform
import io


def formate_abs_path(path: str) -> str:
    """这主要针对windows环境，输入的绝对路径可能不包含盘符，这里进行补充
    主要是用于打印效果
    如果不是windows环境，直接返回path，相当于没有调用这个函数

    Parameters
    ----------
    path : str
        待转换的路径

    Returns
    -------
    str
        增加了盘符的路径
    """
    if platform.system() != "Windows":
        return path
    if not os.path.isabs(path):
        return path
    need_add = len(path) < 3 or path[1] != ":"
    # 处理反斜杠, 保证路径的正确性
    path = path.replace("/", "\\")
    if need_add:
        return os.path.join(os.getcwd()[:2], path)
    return path


def check_dir_and_create(path: str) -> str:
    """检查路径是否是一个绝对路径文件夹
    检查逻辑是：
    1. 检查path输入是否是一个字符串，如果不是，抛出ValueError
    2. 检查path是否是一个绝对路径，如果不是，转换为绝对路径
    3. 创建文件夹，如果文件夹已经存在，则不做任何操作
    4. 判断目标文件夹的权限，如果没有写入权限，则抛出IOError

    Parameters
    ----------
    path : str
        待检查的路径

    Returns
    -------
    str
        将path转换为绝对路径后返回

    Raises
    ------
    ValueError
        path不是一个字符串
    IOError
        path没有写权限
    """
    if not isinstance(path, str):
        raise ValueError("path must be a string")
    if not os.path.isabs(path):
        path = os.path.abspath(path)
    # 如果创建失败，也是抛出IOError
    try:
        os.makedirs(path, exist_ok=True)
    except Exception as e:
        raise IOError(f"create path: {path} failed, error: {e}")
    if not os.access(path, os.W_OK):
        raise IOError(f"no write permission for path: {path}")
    return path


def get_file_hash_path(file_path: str) -> str:
    """计算并返回给定文件的SHA-256哈希值。"""
    import hashlib

    hash_sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:  # 以二进制读取模式打开文件
        while chunk := f.read(8192):  # 读取文件的小块进行处理
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()


def get_file_hash_numpy_array(array) -> str:
    """计算并返回给定NumPy数组的SHA-256哈希值。"""
    import hashlib

    hash_sha256 = hashlib.sha256()
    # 将NumPy数组转换为字节串，然后更新哈希值
    hash_sha256.update(array.tobytes())
    return hash_sha256.hexdigest()


def get_file_hash_pil(image) -> str:
    """计算并返回给定PIL.Image对象的SHA-256哈希值。"""
    import hashlib

    hash_sha256 = hashlib.sha256()
    # 将图像转换为字节数据
    with io.BytesIO() as buffer:
        image.save(buffer, format="PNG")  # 可以选择其他格式，如'JPEG'
        hash_sha256.update(buffer.getvalue())
    return hash_sha256.hexdigest()


def get_text_sha256_hash(text):
    """计算并返回给定文本的SHA-256哈希值。"""
    import hashlib

    # 使用hashlib库创建一个sha256哈希对象
    hash_object = hashlib.sha256()
    # 对输入的文本进行编码，因为hashlib需要的是字节对象
    hash_object.update(text.encode())
    # 获取十六进制格式的哈希值
    hex_dig = hash_object.hexdigest()
    return hex_dig
