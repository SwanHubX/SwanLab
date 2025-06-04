"""
@author: cunyue
@file: models.py
@time: 2025/6/4 15:32
@description: 定义日志备份模型，方便序列化和反序列化操作
"""

import json
import os.path
from typing import Optional, List, Literal, Tuple

from pydantic import BaseModel as PydanticBaseModel
from swankit.callback import ColumnInfo, ColumnConfig, RuntimeInfo, MetricInfo
from swankit.callback.models import ColumnClass, SectionType, YRange
from swankit.core import ChartReference

from swanlab.api.upload import FileModel
from swanlab.api.upload.model import ColumnModel, LogModel
from swanlab.log.type import LogData


class BaseModel(PydanticBaseModel):
    model_type: str = "BaseModel"

    def __getitem__(self, key):
        return getattr(self, key)

    def to_backup(self) -> str:
        """
        将模型转换为 JSON 字符串
        """
        return json.dumps(
            {
                "model_type": self.model_type,
                "data": self.model_dump(exclude_none=True, exclude={"model_type"}),
            }
        )

    @classmethod
    def from_backup(cls, json_str: str):
        """
        从 JSON 字符串创建模型实例
        """
        data = json.loads(json_str)
        model_type = data["model_type"]
        if model_type not in backup_models:
            raise ValueError(f"Unsupported model type: {model_type}")
        return cls.model_validate(
            {
                "model_type": model_type,
                **data["data"],
            }
        )


class Project(BaseModel):
    model_type: Literal["PROJECT"] = "PROJECT"
    name: Optional[str]  # 项目名称
    workspace: Optional[str]  # 工作空间名称
    public: Optional[bool]  # 项目是否公开


class Experiment(BaseModel):
    model_type: Literal["EXPERIMENT"] = "EXPERIMENT"
    name: Optional[str]  # 实验名称
    description: Optional[str]  # 实验描述
    tags: Optional[List[str]]  # 实验标签


class Log(BaseModel):
    model_type: Literal["LOG"] = "LOG"
    create_time: str  # 日志创建时间
    message: str  # 日志消息内容
    epoch: Optional[int]  # 日志所属的训练轮次
    level: Literal['INFO', 'WARN', 'ERROR']  # 日志级别

    @classmethod
    def from_log_data(cls, log_data: LogData) -> List['Log']:
        """
        从 LogData 对象创建 Log 实例
        """
        level = 'INFO' if log_data['type'] == 'stdout' else 'WARN'
        l = []
        for item in log_data['contents']:
            l.append(
                cls.model_validate(
                    {
                        "create_time": item["create_time"],
                        "epoch": item.get("epoch"),
                        "message": item["message"],
                        "level": level,
                    }
                )
            )
        return l

    def to_log_model(self) -> LogModel:
        return LogModel(
            level=self.level,
            contents=[{"create_time": self.create_time, "epoch": self.epoch, "message": self.message}],
        )


class Column(BaseModel):
    model_type: Literal["COLUMN"] = "COLUMN"

    key: str  # 列的唯一标识符
    kid: str  # 列的唯一标识符（可能是一个ID或其他标识）
    name: Optional[str]  # 列的名称
    cls: ColumnClass  # 列的类型，可能是自定义列或系统生成列
    column_type: str  # 列对应的列类型
    chart_reference: ChartReference  # 列对应图表的参考系，可能是步数或时间
    section_name: Optional[str]  # 列所属的组名
    section_type: SectionType  # 列的组类型，可能是自定义组或系统组
    section_sort: Optional[int]  # 列在组中的排序位置
    error: Optional[dict]  # 列的错误信息，如果有的话

    # ---------------------------------- ColumnConfig ----------------------------------
    y_range: YRange  # y轴范围配置
    chart_name: Optional[str]  # 图表名称
    chart_index: Optional[str]  # 图表索引
    metric_name: Optional[str]  # 指标名称
    metric_color: Optional[Tuple[str, str]]  # 指标颜色

    @classmethod
    def from_column_info(cls, column_info: ColumnInfo):
        """
        从 ColumnInfo 对象创建 Column 实例
        """
        return cls.model_validate(
            {
                "key": column_info.key,
                "kid": column_info.kid,
                "name": column_info.name,
                "cls": column_info.cls,
                "column_type": column_info.chart_type.value.column_type,
                "chart_reference": column_info.chart_reference,
                "section_name": column_info.section_name,
                "section_type": column_info.section_type,
                "section_sort": column_info.section_sort,
                "error": column_info.error.dict() if column_info.error else None,
                "y_range": column_info.config.y_range if column_info.config else None,
                "chart_name": column_info.config.chart_name if column_info.config else None,
                "chart_index": column_info.config.chart_index if column_info.config else None,
                "metric_name": column_info.config.metric_name if column_info.config else None,
                "metric_color": column_info.config.metric_color if column_info.config else None,
            }
        )

    def to_column_model(self) -> ColumnModel:
        """
        将 Column 实例转换为 ColumnModel 实例
        """
        return ColumnModel(
            key=self.key,
            name=self.name,
            cls=self.cls,
            typ=self.column_type,
            section_name=self.section_name,
            section_type=self.section_type,
            config=ColumnConfig(
                y_range=self.y_range,
                chart_name=self.chart_name,
                chart_index=self.chart_index,
                metric_name=self.metric_name,
                metric_color=self.metric_color,
            ),
            error=self.error,
        )


class Runtime(BaseModel):
    model_type: Literal["RUNTIME"] = "RUNTIME"
    conda_filename: Optional[str]  # Conda 环境文件名
    requirements_filename: Optional[str]  # Python requirements 文件名
    metadata_filename: Optional[str]  # 系统元数据名
    config_filename: Optional[str]  # 用户自定义配置文件名

    @classmethod
    def from_runtime_info(cls, runtime_info: RuntimeInfo):
        """
        从 RuntimeInfo 对象创建 Runtime 实例
        """
        return cls.model_validate(
            {
                "conda_filename": runtime_info.conda.name if runtime_info.conda else None,
                "requirements_filename": runtime_info.requirements.name if runtime_info.requirements else None,
                "metadata_filename": runtime_info.metadata.name if runtime_info.metadata else None,
                "config_filename": runtime_info.config.name if runtime_info.config else None,
            }
        )

    def to_file_model(self, file_dir) -> FileModel:
        """
        将 Runtime 实例转换为 RuntimeInfo 实例
        """
        return FileModel(
            conda=(
                open(os.path.join(file_dir, self.conda_filename), "r", encoding="utf-8").read()
                if self.conda_filename
                else None
            ),
            requirements=(
                open(os.path.join(file_dir, self.requirements_filename), "r", encoding="utf-8").read()
                if self.requirements_filename
                else None
            ),
            metadata=(
                json.loads(open(os.path.join(file_dir, self.metadata_filename), "r", encoding="utf-8").read())
                if self.metadata_filename
                else None
            ),
            config=(
                json.loads(open(os.path.join(file_dir, self.conda_filename), "r", encoding="utf-8").read())
                if self.conda_filename
                else None
            ),
        )


class Metric(BaseModel):
    model_type: Literal['METRIC'] = "METRIC"
    name: str  # 指标名称
    value: float  # 指标值
    step: Optional[int]  # 指标对应的步数或时间戳

    @classmethod
    def from_metric_info(cls, metric_info: MetricInfo):
        pass


backup_models = {
    "PROJECT": Project,
    "EXPERIMENT": Experiment,
    "LOG": Log,
    "COLUMN": Column,
    "RUNTIME": Runtime,
    "METRIC": Metric,
}
