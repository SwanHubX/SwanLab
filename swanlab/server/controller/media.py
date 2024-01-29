#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-29 14:12:55
@File: swanlab\server\controller\media.py
@IDE: vscode
@Description:
    媒体数据解析
    - 音频
    - ...
"""

import os
import wave
import numpy as np

# 响应相关
from ..module.resp import (
    SUCCESS_200,
    DATA_ERROR_500,
    CONFLICT_409,
    NOT_FOUND_404,
)

# 路径相关
from ..settings import (
    SWANLOG_DIR,
)


# ---------------------------------- 工具函数 ----------------------------------


def __get_wav_header_info(file_path):
    """获取 wav 文件的头部信息

    Parameters
    ----------
    file_path : str
        wav 文件路径

    Returns
    -------
    dict
        - chunk_id : 4个字节,代表文件标识符,通常为 "RIFF"
        - chunk_size : 4个字节,表示文件大小,包括头部以及音频数据部分的大小
        - format : 4个字节,代表文件格式标识符,通常为 "WAVE"
        - sub_chunk_1_id : 4个字节,代表子块1标识符,通常为 "fmt "
        - sub_chunk_1_size : 4个字节,表示子块1的大小,通常为16(表示PCM格式)
        - audio_format : 2个字节,表示音频格式,常见的值有1(表示PCM),其他值表示压缩格式
        - num_channels : 2个字节,表示声道数
        - sample_rate : 4个字节,表示采样率,常见值44100、48000
        - byte_rate : 4个字节,表示字节率,即每秒钟数据传输速率
        - block_align : 2个字节,表示块对齐,即每个采样所占的字节数,计算方法为NumChannels * BitsPerSample / 8
        - bits_per_sample : 2个字节,表示每个样本的位深度,即用来表示每个采样值的二进制位数,常见用8、16、24
    """

    header_info = {}
    with open(file_path, "rb") as f:
        # 读取头部信息
        riff_chunk = f.read(12)
        fmt_chunk = f.read(24)

        # 解析RIFF chunk
        riff_chunk_id = riff_chunk[:4]
        riff_chunk_size = int.from_bytes(riff_chunk[4:8], byteorder="little")
        riff_format = riff_chunk[8:12]

        # 解析fmt chunk
        fmt_chunk_id = fmt_chunk[:4]
        fmt_chunk_size = int.from_bytes(fmt_chunk[4:8], byteorder="little")
        audio_format = int.from_bytes(fmt_chunk[8:10], byteorder="little")
        num_channels = int.from_bytes(fmt_chunk[10:12], byteorder="little")
        sample_rate = int.from_bytes(fmt_chunk[12:16], byteorder="little")
        byte_rate = int.from_bytes(fmt_chunk[16:20], byteorder="little")
        block_align = int.from_bytes(fmt_chunk[20:22], byteorder="little")
        bits_per_sample = int.from_bytes(fmt_chunk[22:24], byteorder="little")

        # 存储头部信息到字典中
        header_info = {
            "chunk_id": riff_chunk_id.decode("utf-8"),
            "chunk_size": riff_chunk_size,
            "format": riff_format.decode("utf-8"),
            "sub_chunk_1_id": fmt_chunk_id.decode("utf-8"),
            "sub_chunk_1_size": fmt_chunk_size,
            "audio_format": audio_format,
            "num_channels": num_channels,
            "sample_rate": sample_rate,
            "byte_rate": byte_rate,
            "block_align": block_align,
            "bits_per_sample": bits_per_sample,
        }

    return header_info


def __get_wav_duration(path: str, sample_rate: float):
    """获取 wav 音频的时长

    Parameters
    ----------
    path : str
        wav 文件路径
    sample_rate : float
        wav 的采样速率

    Returns
    -------
    duration : float
        音频时长，单位秒
    """

    # 打开.wav文件
    with wave.open(path, "r") as wav_file:
        # 获取帧数
        frames = wav_file.getnframes()
        # 获取采样率
        sample_rate = wav_file.getframerate()
        # 计算音频时长（以秒为单位）
        duration = frames / float(sample_rate)

    return duration


# ---------------------------------- 音频相关 ----------------------------------


MAX_LENGTH_PER_SECOND = 200


def get_audio_data(path: str):
    """解析并获取音频数据

    Parameters
    ----------
    path : str
        音频文件相对于 swanlog 目录的相对路径

    Returns
    -------
    data : dict
        - path : string
            音频文件路径
        - header_info : dict
            音频头部信息
        - data : list
            音频数据
        - duration : float
            音频时长
    """

    # TODO 处理路径的拼接
    path = os.path.join(SWANLOG_DIR, path)
    # 获取头部信息
    header_info = __get_wav_header_info(path)
    # 获取音频时长
    duration = round(__get_wav_duration(path, header_info["sample_rate"]), 2)
    # 获取音频数据
    with open(path, "rb") as f:
        # 跳过WAV文件头部
        f.read(44)
        # 读取音频数据
        data = np.fromfile(f, dtype=np.int16)

    # 对data进行取样
    max_length = int(duration * MAX_LENGTH_PER_SECOND)
    if len(data) > max_length:
        sampling_ratio = len(data) // max_length  # 计算取样比例
        data = data[::sampling_ratio]  # 进行取样

    return SUCCESS_200(
        {
            "path": path,
            "header_info": header_info,
            "data": data.tolist(),
            "duration": duration,
        }
    )
