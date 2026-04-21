"""
@author: cunyue
@file: builder.py
@time: 2026/4/21 16:56
@description: 构建记录，核心在于一些字段约定
"""

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord

__all__ = ["build_finish_record", "build_start_record"]

START_RECORD_NUM = -1
"""
约定启动记录的num为-1
"""
FINISH_RECORD_NUM = -2
"""
约定结束记录的num为-2
"""


def build_start_record(start_record: StartRecord):
    """
    构建启动记录
    """
    return Record(num=START_RECORD_NUM, start=start_record)


def build_finish_record(finish_record: FinishRecord):
    """
    构建结束记录
    """
    return Record(num=FINISH_RECORD_NUM, finish=finish_record)
