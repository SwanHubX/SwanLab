"""
@author: cunyue
@file: _writer.py
@time: 2026/3/14
@description: Config 文件写入层
将内存中的 config 序列化为后端约定的 YAML 格式后，通过 safe_write 落盘。
"""

from pathlib import Path

import yaml

from swanlab.sdk.internal.pkg.fs import safe_write

__all__ = ["write_config"]


def write_config(path: Path, config: dict, sort_map: dict) -> None:
    """
    将 config 序列化为 {key: {value, desc, sort}} 格式并写入 YAML 文件。

    每次调用均全量覆写（INIT 和 PATCH 均如此），消费方以最新文件内容为准。

    :param path:     目标文件路径（config.yaml）
    :param config:   内部存储的原始 {key: value} dict（value 已经过 parse()）
    :param sort_map: key → sort index 映射，用于还原插入顺序
    """
    formatted = {key: {"value": value, "desc": "", "sort": sort_map.get(key, 0)} for key, value in config.items()}
    content = yaml.safe_dump(formatted, allow_unicode=True, default_flow_style=False)
    safe_write(path, content)
