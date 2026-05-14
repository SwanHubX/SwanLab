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
from swanlab.proto.swanlab.save.v1.save_pb2 import SavePolicy
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogLevel

from .bimap import BiMap

__all__ = ["resume", "medium", "state", "level", "policy"]


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
        "aborted": RunState.RUN_STATE_ABORTED,
    }
)
"""RunState 枚举适配器"""


column = BiMap(
    {
        "FLOAT": ColumnType.COLUMN_TYPE_SCALAR,
        "IMAGE": ColumnType.COLUMN_TYPE_IMAGE,
        "AUDIO": ColumnType.COLUMN_TYPE_AUDIO,
        "VIDEO": ColumnType.COLUMN_TYPE_VIDEO,
        "TEXT": ColumnType.COLUMN_TYPE_TEXT,
        "ECHARTS": ColumnType.COLUMN_TYPE_ECHARTS,
        "OBJECT3D": ColumnType.COLUMN_TYPE_OBJECT3D,
        "MOLECULE": ColumnType.COLUMN_TYPE_MOLECULE,
    }
)
"""ColumnType 枚举适配器，映射本地proto枚举与云端类型枚举"""


medium = BiMap(
    {
        "text": ColumnType.COLUMN_TYPE_TEXT,
        "image": ColumnType.COLUMN_TYPE_IMAGE,
        "audio": ColumnType.COLUMN_TYPE_AUDIO,
        "video": ColumnType.COLUMN_TYPE_VIDEO,
        "echarts": ColumnType.COLUMN_TYPE_ECHARTS,
        "object3d": ColumnType.COLUMN_TYPE_OBJECT3D,
        "molecule": ColumnType.COLUMN_TYPE_MOLECULE,
    }
)
"""媒体类型名称与 ColumnType 枚举的双向映射，同时作为存储目录名映射"""


level = BiMap(
    {
        "INFO": LogLevel.LOG_LEVEL_INFO,
        "ERROR": LogLevel.LOG_LEVEL_ERROR,
    }
)
"""LogLevel 枚举适配器，映射本地proto枚举与云端日志格式"""

policy = BiMap(
    {
        "now": SavePolicy.SAVE_POLICY_NOW,
        "end": SavePolicy.SAVE_POLICY_END,
        "live": SavePolicy.SAVE_POLICY_LIVE,
    }
)

"""Save Policy 文件保存策略适配器，映射 save 时文件上报策略"""
