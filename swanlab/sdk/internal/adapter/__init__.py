"""
@author: cunyue
@file: __init__.py
@time: 2026/3/15 00:56
@description: 适配器，将一个类的接口转换成客户希望的另外一个接口。
举个例子，用户传入的resume参数是"must"、"allow"、"never"，但我们内部使用Protobuf枚举
两者需要适配和转换
"""

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.run.v1.run_pb2 import ResumeMode, RunState

from .bimap import BiMap

__all__ = ["resume", "column_type", "state"]

resume = BiMap(
    {
        "allow": ResumeMode.RESUME_MODE_ALLOW,
        "must": ResumeMode.RESUME_MODE_MUST,
        "never": ResumeMode.RESUME_MODE_NEVER,
    }
)
"""ResumeMode 枚举适配器"""

state = BiMap(
    {
        "success": RunState.RUN_STATE_FINISHED,
        "crashed": RunState.RUN_STATE_CRASHED,
        "aborted": RunState.RUN_STATE_STOPPED,
    }
)
"""RunState 枚举适配器"""

column_type = BiMap(
    {
        "scalar": ColumnType.COLUMN_TYPE_FLOAT,
        "text": ColumnType.COLUMN_TYPE_TEXT,
        "image": ColumnType.COLUMN_TYPE_IMAGE,
        "audio": ColumnType.COLUMN_TYPE_AUDIO,
        "video": ColumnType.COLUMN_TYPE_VIDEO,
        "echarts": ColumnType.COLUMN_TYPE_ECHARTS,
    }
)
"""ColumnType 枚举适配器，同时也会作为路径映射"""
