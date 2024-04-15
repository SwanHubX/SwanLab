#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-17 15:13:31
@File: swanlab/server/controller/utils/tag.py
@IDE: vscode
@Description:
    tag相关处理函数
"""
from typing import List
import json
import ujson
import math
import os

# tag 总结文件名
TAG_SUMMARY_FILE = "_summary.json"
# logs 目录下的配置文件
LOGS_CONFIGS = [TAG_SUMMARY_FILE]


def read_tag_data(file_path: str) -> List[dict]:
    """读取某一个tag文件数据，数据依据其中的index字段排序

    Parameters
    ----------
    file_path : str
        tag文件路径，绝对路径
    """
    # 如果文件内容为空，返回空列表
    if os.path.getsize(file_path) == 0:
        return []
    # 读取文件内容，文件内部本质上是字符串一堆json格式的数据，每一行是一个json并且有换行符分隔，需要拿到然后解析为list<dict>
    data = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        # 解析为list<dict>
        for i in range(len(lines)):
            if len(lines[i]):
                data.append(json.loads(lines[i]))
        return data


def get_tag_files(tag_path: str, exclude=None) -> List[str]:
    """
    获取实验数据，并且做向下兼容，在v0.2.4版本以前的实验日志数据将转换为新的格式

    Parameters
    ----------
    tag_path : str
        tag的路径
    exclude : List[str]
        需要排除的文件列表

    Returns
    -------
    List[str]
        tag文件列表
    """
    # 降序排列，最新的数据在最前面
    if exclude is None:
        exclude = []
    files: list = os.listdir(tag_path)
    previous_logs = [f for f in files if f.endswith(".json") and f not in exclude]
    current_logs = [f for f in files if f.endswith(".log")]
    current_logs.sort()
    # COMPAT 如果目标文件夹不存在*.log文件但存在*.json文件，说明是之前的日志格式，需要转换
    if len(current_logs) == 0 and len(previous_logs) > 0:
        for file in previous_logs:
            with open(os.path.join(tag_path, file), "r") as f:
                data = ujson.load(f)
                with open(os.path.join(tag_path, file.replace(".json", ".log")), "a") as i:
                    for d in data["data"]:
                        i.write(ujson.dumps(d) + "\n")
        current_logs = [f for f in os.listdir(tag_path) if f.endswith(".log")]
        current_logs.sort()
    return current_logs


def _sample_a_bucket(bucket: List[dict]) -> dict:
    """
    对一个bucket进行采样，返回采样后的数据
    :param bucket: 一个桶中的信息
    :return: 采样后的数据
    """
    if len(bucket) <= 2:
        return bucket[0]
    # 计算bucket头尾直线表达式
    tail = bucket[0]
    head = bucket[-1]
    # k，b
    k = (head["data"] - tail["data"]) / (head["index"] - tail["index"])
    b = head["data"] - k * head["index"]
    # 计算每个点到直线的距离
    max_distance_p = (0, 0)  # (d, i)
    for i in range(1, len(bucket) - 1):
        distance = abs(k * bucket[i]["index"] + b - bucket[i]["data"]) / math.sqrt(k ** 2 + 1)
        if distance > max_distance_p[0]:
            max_distance_p = (distance, i)
    return bucket[max_distance_p[1]]


def _calculate_bucket_capacity(total_data_points, target_bucket_count) -> List[int]:
    """
    计算每个桶的容量
    :param total_data_points: 总数据点数
    :param target_bucket_count: 目标桶数量
    :return: 每个桶的容量
    """
    # 计算每个桶的容量
    bucket_capacity = total_data_points // target_bucket_count  # 整除部分
    remaining_points = total_data_points % target_bucket_count  # 余数

    # 将余数均匀地分配到前几个桶中
    bucket_capacity_list = [bucket_capacity] * target_bucket_count
    for i in range(remaining_points):
        bucket_capacity_list[i] += 1

    return bucket_capacity_list


def lttb(data: List[dict], threshold: int = 1500):
    """
    LTTB算法，对数据进行降采样
    最后保留头尾数据，因此桶的数量应该是实际数量-2
    桶大小向下取整，保留第一个和最后一个点
    降采样将优先发生在尾部
    :param data: 数据
    :param threshold: 采样后的数据长度，也是降采样阈值
    """
    if len(data) <= threshold:
        return data
    bucket_n = threshold - 2
    _data = data[1:-1]
    # 需要保留第一、最后一个点
    sampled = [data[0]]
    # 计算每一个桶的容量，每个桶容量不尽相同
    buckets_capacity = _calculate_bucket_capacity(len(_data), bucket_n)
    # 当前使用的bucket
    now_bucket = []
    # 当前bucket的数量
    now_bucket_n = 0
    for _d in _data:
        now_bucket.append(_d)
        if len(now_bucket) < buckets_capacity[now_bucket_n]:
            continue
        sampled.append(_sample_a_bucket(now_bucket))
        now_bucket = []
        now_bucket_n += 1
        if now_bucket_n == bucket_n - 1:
            break
    sampled.append(data[-1])
    return sampled
