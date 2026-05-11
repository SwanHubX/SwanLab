"""SwanLab Core settings."""

import os

from pydantic import BaseModel, ConfigDict, Field


def section_rule_index_factory() -> int:
    # 向下兼容旧版本环境变量
    return int(os.environ.get("SWANLAB_SECTION_RULE_IDX", "0"))


class SaveSettings(BaseModel):
    """swanlab.save() 文件上传限制。"""

    multipart_threshold: int = Field(default=100 * 1024 * 1024)
    """超过此大小使用分片上传（默认 100 MB）。"""

    part_size: int = Field(default=32 * 1024 * 1024)
    """分片上传时每片大小（默认 32 MB）。"""

    max_file_size: int = Field(default=50 * 1024 * 1024 * 1024)
    """单文件上传大小上限（默认 50 GB）。"""

    max_batch_size: int = Field(default=100)
    """每批上传文件数量上限（默认 100）。"""

    max_total_size: int = Field(default=50 * 1024 * 1024 * 1024)
    """总文件大小上限（默认 50 GB）。"""

    model_config = ConfigDict(frozen=True)


class MediaSettings(BaseModel):
    """媒体记录（Image/Audio/Video 等）限制。"""

    pass


class RecordSettings(BaseModel):
    """Record 上传批次限制。"""

    batch_interval: float = Field(default=5.0, gt=0)
    """
    Batch interval (seconds) for the Transport upload thread.
    Default 5.0s. Use a smaller value (e.g. 0.5) in Converter scenarios for higher throughput.
    """

    batch_size: int = Field(default=10_000, gt=0)
    """
    Batch size of records per HTTP request in Dispatch. Default 10000. 
    """


class CoreSettings(BaseModel):
    """
    Slash index used to split a metric key into section name and metric name.
    0 means the first slash, 1 means the second slash, -1 means the last slash.
    Out-of-range values wrap via ``idx % slash_count``.
    The split position is fixed when settings are loaded for each run.
    """

    section_rule_index: int = Field(default_factory=section_rule_index_factory)

    model_config = ConfigDict(frozen=True)

    Save: type = SaveSettings
    Media: type = MediaSettings
    Record: type = RecordSettings

    save: SaveSettings = Field(default_factory=SaveSettings)
    media: MediaSettings = Field(default_factory=MediaSettings)
    record: RecordSettings = Field(default_factory=RecordSettings)
