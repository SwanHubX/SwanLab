"""SwanLab Core settings."""

from pydantic import BaseModel, ConfigDict, Field

from .compat import section_rule_index_factory


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
    record_interval: float = Field(default=1.5, gt=0)
    """
    Batch interval (seconds) for the Transport upload thread.
    Default 1.5s. Use a smaller value (e.g. 0.5) in Converter scenarios for higher throughput.
    """
    media_max_size: int = Field(default=50 * 1024**2, ge=10 * 1024**2, le=500 * 1024**2)
    """
    Maximum size in bytes of a single media file. Default 50 MiB, range 10 MiB to 500 MiB.

    Applies uniformly to every media type (image, audio, video, ...); records whose items
    exceed the limit are dropped with a warning. In ``skip_store`` mode the limit also caps
    the inlined payload of a single media record item.
    """
    save_split: int = Field(default=100 * 1024 * 1024)
    """
    File size threshold in bytes for multipart upload. Default 100 MiB.
    """
    save_size: int = Field(default=50 * 1024 * 1024 * 1024)
    """
    Maximum saved size per file. Default 50 GiB.
    """
    save_part: int = Field(default=32 * 1024 * 1024, gt=0)
    """
    Multipart upload part size in bytes. Default 32 MiB.
    """
    save_batch: int = Field(default=100)
    """
    Maximum number of files per save upload batch. Default 100.
    """
    skip_store: bool = Field(default=False)
    """
    Skip all SDK-owned local files and upload to the cloud only (requires ``online`` mode).

    Media contents and internal texts (config/metadata/requirements/conda) are inlined
    into records and uploaded directly; ``run.save()`` files are read from their original
    path without local links.
    ``swanlab sync`` / ``swanlab watch`` do not apply; a crashed run cannot be recovered.
    """
    model_config = ConfigDict(frozen=True)
