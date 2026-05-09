"""SwanLab Core settings."""

import os

from pydantic import BaseModel, ConfigDict, Field


def section_rule_index_factory() -> int:
    # 向下兼容旧版本环境变量
    return int(os.environ.get("SWANLAB_SECTION_RULE_IDX", "0"))


class CoreSettings(BaseModel):
    section_rule_index: int = Field(default_factory=section_rule_index_factory)
    """
    Slash index used to split a metric key into section name and metric name.
    0 means the first slash, 1 means the second slash, -1 means the last slash.
    Out-of-range values wrap via ``idx % slash_count``.
    The split position is fixed when settings are loaded for each run.
    """

    batch_interval: float = Field(default=5.0)
    """
    Batch interval (seconds) for the Transport upload thread.
    Default 5.0s. Use a smaller value (e.g. 0.5) in Converter scenarios for higher throughput.
    """

    max_records_per_request: int = Field(default=10_000)
    """
    Maximum number of records per HTTP request in Dispatch. Default 10000. 
    """

    finish_join_timeout: int = Field(default=30)
    """
    Timeout (seconds) when waiting for the Transport thread to exit on finish. Default 30.
    """

    model_config = ConfigDict(frozen=True)
