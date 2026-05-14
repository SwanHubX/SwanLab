"""SwanLab Core settings."""

import os

from pydantic import BaseModel, ConfigDict, Field


def section_rule_index_factory() -> int:
    # 向下兼容旧版本环境变量
    return int(os.environ.get("SWANLAB_SECTION_RULE_IDX", "0"))


class CoreSettings(BaseModel):
    """
    Slash index used to split a metric key into section name and metric name.
    0 means the first slash, 1 means the second slash, -1 means the last slash.
    Out-of-range values wrap via ``idx % slash_count``.
    The split position is fixed when settings are loaded for each run.
    """

    section_rule: int = Field(default_factory=section_rule_index_factory)
    """
    Section rule for auto generated section name.
    """
    record_batch: int = Field(default=10_000, gt=0, lt=100_000)
    """
    Batch size of records per HTTP request in Dispatch. Default 10000. 
    """
    record_interval: float = Field(default=5.0, gt=0)
    """
    Batch interval (seconds) for the Transport upload thread.
    Default 5.0s. Use a smaller value (e.g. 0.5) in Converter scenarios for higher throughput.
    """
    save_split: int = Field(default=100 * 1024 * 1024)
    """
    File size threshold in bytes for multipart upload. Default 100 MiB.
    """
    save_size: int = Field(default=50 * 1024 * 1024 * 1024)
    """
    Maximum saved size per file. Default 50 GiB.
    """
    save_part: int = Field(default=32 * 1024 * 1024)
    """
    Multipart upload part size in bytes. Default 32 MiB.
    """
    save_batch: int = Field(default=100)
    """
    Maximum number of files per save upload batch. Default 100.
    """
    model_config = ConfigDict(frozen=True)
