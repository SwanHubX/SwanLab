"""
@author: cunyue
@file: _writer.py
@time: 2026/3/14
@description: Config 文件写入层
将内存中的 config 序列化为后端约定的 YAML 格式后，通过 safe_write 落盘。
"""

from pathlib import Path

import yaml

from swanlab.sdk.internal.pkg import fs

__all__ = ["dump_config", "format_config", "write_config"]


def format_config(config: dict, sort_map: dict) -> dict:
    """
    将内部存储的 {key: value} 归一化为后端约定的 {key: {value, desc, sort}} 结构。

    :param config:  内部存储的原始 {key: value} dict（value 已经过 parse()）
    :param sort_map: key → sort index 映射，用于还原插入顺序
    """
    return {key: {"value": value, "desc": "", "sort": sort_map.get(key, 0)} for key, value in config.items()}


def dump_config(content: dict) -> str:
    """
    将归一化后的 config 结构序列化为 YAML 文本。
    落盘（write_config）与 skip_store 下的内联上传共用同一份编码，避免两种模式云端 config 结构漂移。

    :param content: format_config 产出的归一化结构
    """
    return yaml.safe_dump(content, allow_unicode=True, default_flow_style=False)


def write_config(path: Path, config: dict, sort_map: dict) -> None:
    """
    将 config 序列化为 {key: {value, desc, sort}} 格式并写入 YAML 文件。

    每次调用均全量覆写（INIT 和 PATCH 均如此），消费方以最新文件内容为准。

    :param path:     目标文件路径（config.yaml）
    :param config:   内部存储的原始 {key: value} dict（value 已经过 parse()）
    :param sort_map: key → sort index 映射，用于还原插入顺序
    """
    fs.safe_write(path, dump_config(format_config(config, sort_map)))
