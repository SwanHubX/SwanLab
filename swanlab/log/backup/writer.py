"""
@author: cunyue
@file: writer.py
@time: 2025/6/4 17:49
@description: 写入器，因为在备份写入本地的时候，local 模式写入的内容可能与备份重合，因此将这部分逻辑封装出来
"""

import os

from swanlab.toolkit import RuntimeInfo, MetricInfo


def write_runtime_info(files_dir: str, runtime_info: RuntimeInfo):
    """
    写入运行时信息到文件
    :param files_dir: 文件目录
    :param runtime_info: 运行时信息
    """
    # 更新运行时信息
    if runtime_info.requirements is not None:
        runtime_info.requirements.write(files_dir)
    if runtime_info.metadata is not None:
        runtime_info.metadata.write(files_dir)
    if runtime_info.config is not None:
        runtime_info.config.write(files_dir)
    if runtime_info.conda is not None:
        runtime_info.conda.write(files_dir)


def write_media_buffer(metric_info: MetricInfo):
    """
    写入媒体缓冲区到指定目录
    """
    if metric_info.metric_buffers is None:
        return
    for i, r in enumerate(metric_info.metric_buffers):
        if r is None:
            continue
        # 组合路径
        path = os.path.join(metric_info.swanlab_media_dir, metric_info.column_info.kid)
        os.makedirs(path, exist_ok=True)
        # 写入数据
        with open(os.path.join(path, metric_info.metric["data"][i]), "wb") as f:
            f.write(r.getvalue())
