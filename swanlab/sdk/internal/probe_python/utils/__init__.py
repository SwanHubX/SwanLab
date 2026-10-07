"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/07
@description: probe_python 工具函数
"""

from pathlib import Path

from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.pkg import fs


def make_save_record(name: str, save_type: SaveType, content: str, target: Path, skip_store: bool) -> SaveRecord:
    """构建 SaveRecord。

    当 skip_store 为 True 时，content 写入 payload 且不落盘；否则写入目标文件并引用路径。
    """
    if skip_store:
        return SaveRecord(name=name, type=save_type, payload=content.encode("utf-8"))
    fs.safe_write(target, content)
    return SaveRecord(name=name, type=save_type, source_path=target.absolute().as_posix())
