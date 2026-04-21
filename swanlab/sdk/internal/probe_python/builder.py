"""
@author: cunyue
@file: builder.py
@time: 2026/4/21 15:40
@description: 构建硬件监控记录
"""

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.system.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord

METADATA_RECORD_NUM = -4
"""
约定元信息记录的num为-4
"""
REQUIREMENTS_RECORD_NUM = -5
"""
约定依赖记录的num为-5
"""
CONDA_RECORD_NUM = -6
"""
约定conda环境记录的num为-6
"""


def build_metadata_record(ts: Timestamp):
    """
    构建元信息记录
    """
    return Record(num=METADATA_RECORD_NUM, metadata=MetadataRecord(timestamp=ts))


def build_requirements_record(ts: Timestamp):
    """
    构建依赖记录
    """
    return Record(num=REQUIREMENTS_RECORD_NUM, requirements=RequirementsRecord(timestamp=ts))


def build_conda_record(ts: Timestamp):
    """
    构建conda环境记录
    """
    return Record(num=CONDA_RECORD_NUM, conda=CondaRecord(timestamp=ts))
