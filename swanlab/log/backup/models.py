"""
@author: cunyue
@file: models.py
@time: 2025/6/4 15:32
@description: 定义日志备份模型，方便序列化和反序列化操作
"""

import json
import os.path
from typing import Optional, List, Literal, Tuple, Union

import yaml
from pydantic import BaseModel as PydanticBaseModel
from swankit.callback import ColumnInfo, ColumnConfig, RuntimeInfo, MetricInfo
from swankit.callback.models import ColumnClass, SectionType, YRange
from swankit.core import ChartReference, MediaBuffer

from swanlab.api.upload import FileModel, ScalarModel
from swanlab.api.upload.model import ColumnModel, LogModel, MediaModel
from swanlab.log.type import LogData


class BaseModel(PydanticBaseModel):

    def __getitem__(self, key):
        return getattr(self, key)

    def to_record(self) -> str:
        """
        将模型转换为 JSON 字符串
        """
        return (
            json.dumps(
                {
                    "model_type": type(self).__name__,
                    "data": self.model_dump(),
                },
                ensure_ascii=False,
            )
            + "\n"
        )

    @classmethod
    def from_record(cls, data: str) -> Tuple[str, 'BaseModel']:
        """
        从 JSON 字符串创建模型实例
        """
        # 最后一个字符串为换行符，需要删除
        data = json.loads(data[:-1])
        model_type = data["model_type"]
        if model_type not in backup_models:
            raise ValueError(f"Unsupported model type: {model_type}")
        return backup_models[model_type].model_validate(data["data"])


class Header(BaseModel):
    backup_type: Literal["DEFAULT"]
    create_time: str  # 备份文件创建时间


class Footer(BaseModel):
    success: bool
    create_time: str


class Project(BaseModel):
    name: Optional[str]  # 项目名称
    workspace: Optional[str]  # 工作空间名称
    public: Optional[bool]  # 项目是否公开


class Experiment(BaseModel):
    name: Optional[str]  # 实验名称
    description: Optional[str]  # 实验描述
    tags: Optional[List[str]]  # 实验标签


class Log(BaseModel):
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


class Runtime(BaseModel):
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
        conda, requirements, metadata, config = None, None, None, None
        if self.conda_filename:
            conda_path = os.path.join(file_dir, self.conda_filename)
            if not os.path.exists(os.path.join(file_dir, self.conda_filename)):
                raise FileNotFoundError(f"Conda file {self.conda_filename} not found in {file_dir}")
            conda = open(conda_path, "r").read()
        if self.requirements_filename:
            requirements_path = os.path.join(file_dir, self.requirements_filename)
            if not os.path.exists(requirements_path):
                raise FileNotFoundError(f"Requirements file {self.requirements_filename} not found in {file_dir}")
            requirements = open(requirements_path, "r").read()
        if self.metadata_filename:
            metadata_path = os.path.join(file_dir, self.metadata_filename)
            if not os.path.exists(metadata_path):
                raise FileNotFoundError(f"Metadata file {self.metadata_filename} not found in {file_dir}")
            metadata = json.loads(open(metadata_path, "r").read())
        if self.config_filename:
            config_path = os.path.join(file_dir, self.config_filename)
            if not os.path.exists(config_path):
                raise FileNotFoundError(f"Config file {self.config_filename} not found in {file_dir}")
            config = yaml.safe_load(open(config_path, "r").read())

        return FileModel(conda=conda, requirements=requirements, metadata=metadata, config=config)


class Column(BaseModel):

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


class Metric(BaseModel):

    @classmethod
    def from_metric_info(cls, metric_info: MetricInfo) -> Union['Scalar', 'Media']:
        """
        指标类型比较特殊，因为上传到后端时分不同接口上传
        此函数根据 MetricInfo 的类型返回对应的 Metric 子类实例
        :param metric_info: MetricInfo 对象
        """
        # 标量类型
        if metric_info.column_info.chart_type == metric_info.column_info.chart_type.LINE:
            return Scalar.model_validate(
                {
                    "metric": metric_info.metric,
                    "key": metric_info.column_info.key,
                    "step": metric_info.metric_step,
                    "epoch": metric_info.metric_epoch,
                }
            )
        buffers_name = []
        if metric_info.metric_buffers is not None:
            for name in metric_info.metric['data']:
                buffers_name.append(name)

        # 媒体类型
        return Media.model_validate(
            {
                "metric": metric_info.metric,
                "key": metric_info.column_info.key,
                "kid": metric_info.column_info.kid,
                "key_encoded": metric_info.column_info.key_encode,
                "step": metric_info.metric_step,
                "epoch": metric_info.metric_epoch,
                "buffers_name": buffers_name if len(buffers_name) else None,
            }
        )


class Scalar(Metric):
    metric: dict  # 标量指标数据，通常是一个字典，包含指标名称和对应的值
    key: str  # 标量指标的唯一标识符
    step: int  # 标量指标的步数
    epoch: int  # 标量指标的训练轮次

    def to_scalar_model(self) -> ScalarModel:
        """
        将 Scalar 实例转换为 ScalarModel 实例
        """
        return ScalarModel(
            metric=self.metric,
            key=self.key,
            step=self.step,
            epoch=self.epoch,
        )


class Media(Metric):
    metric: dict  # 媒体指标数据，通常是一个字典，包含媒体类型和对应的文件路径或URL
    key: str  # 媒体指标的唯一标识符
    kid: int  # 当前实验下，列的唯一id，与保存路径等信息有关，与云端请求无关
    key_encoded: str  # 编码后的键值
    step: int  # 媒体指标的步数
    epoch: int  # 媒体指标的训练轮次
    buffers_name: Optional[List[str]]  # 媒体文件的名称

    def to_media_model(self, media_dir: str) -> MediaModel:
        """
        将 Media 实例转换为 MediaModel 实例
        """
        buffers = []
        # 回复原本的 MediaBuffer 对象
        if self.buffers_name:
            for i, buffer_name in enumerate(self.buffers_name):
                buffer = MediaBuffer()
                buffer.write(open(os.path.join(media_dir, str(self.kid), buffer_name), "rb").read())
                buffer.file_name = "{}/{}".format(self.key_encoded, self.metric["data"][i])
                buffers.append(buffer)

        return MediaModel(
            metric=self.metric,
            key=self.key,
            key_encoded=self.key_encoded,
            step=self.step,
            epoch=self.epoch,
            buffers=buffers if len(buffers) else None,
        )


backup_models = {
    model.__name__: model
    for model in [
        Header,
        Project,
        Experiment,
        Log,
        Runtime,
        Column,
        Scalar,
        Media,
        Footer,
    ]
}


class ModelsParser:
    def __init__(self):
        self._header: Optional[Header] = None
        self._project: Optional[Project] = None
        self._experiment: Optional[Experiment] = None
        self._logs: List[Log] = []
        self._runtime: Runtime = Runtime(
            conda_filename=None,
            requirements_filename=None,
            metadata_filename=None,
            config_filename=None,
        )
        self._columns: List[Column] = []
        self._scalars: List[Scalar] = []
        self._medias: List[Media] = []
        self._footer: Optional[Footer] = None
        self._parsed = False

    def parse_record(self, data: str):
        assert self._parsed, "Must parse records in a context manager"
        record = BaseModel.from_record(data)
        if isinstance(record, Header):
            assert self._header is None, "Header already parsed"
            self._header = record
            return
        if isinstance(record, Project):
            assert self._project is None, "Project already parsed"
            self._project = record
            return
        if isinstance(record, Experiment):
            assert self._experiment is None, "Experiment already parsed"
            self._experiment = record
            return
        if isinstance(record, Log):
            self._logs.append(record)
            return
        if isinstance(record, Runtime):
            if record.conda_filename is not None:
                self._runtime.conda_filename = record.conda_filename
            if record.requirements_filename is not None:
                self._runtime.requirements_filename = record.requirements_filename
            if record.metadata_filename is not None:
                self._runtime.metadata_filename = record.metadata_filename
            if record.config_filename is not None:
                self._runtime.config_filename = record.config_filename
            return
        if isinstance(record, Column):
            self._columns.append(record)
            return
        if isinstance(record, Scalar):
            self._scalars.append(record)
            return
        if isinstance(record, Media):
            self._medias.append(record)
            return
        if isinstance(record, Footer):
            assert self._footer is None, "Footer already parsed"
            self._footer = record
            return

        raise TypeError("Unsupported record type: {}".format(type(record).__name__))

    def __enter__(self):
        assert not self._parsed, "Already parsed records before entering context"
        self._parsed = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # 出现错误则抛出异常
        if exc_type is not None:
            raise exc_val
        # 检查是否所有必要的记录都已解析
        assert self._header is not None, "Header not parsed"
        assert self._project is not None, "Project not parsed"
        assert self._experiment is not None, "Experiment not parsed"

    def get_parsed(
        self,
    ) -> Tuple[
        Header,
        Project,
        Experiment,
        List[Log],
        Runtime,
        List[Column],
        List[Scalar],
        List[Media],
        Optional[Footer],
    ]:
        """
        获取已解析的记录
        """
        assert self._parsed, "Must parse records before getting parsed data"
        # 解析运行时信息，

        return (
            self._header,
            self._project,
            self._experiment,
            self._logs,
            self._runtime,
            self._columns,
            self._scalars,
            self._medias,
            self._footer,
        )
