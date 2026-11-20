"""SwanLab Core settings."""

from typing import Optional

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
    upload_budget_records: int = Field(default=0)
    """
    Upload pipeline budget in records; ``<=0`` means unlimited.
    """
    upload_budget_bytes: int = Field(default=0)
    """
    Upload pipeline budget in bytes; ``<=0`` means unlimited.
    """
    health_probe_interval: float = Field(default=5.0, ge=0)
    """
    Site health probe interval in seconds. ``0`` lets core use its default (5s).
    """
    health_probe_timeout: float = Field(default=2.0, ge=0)
    """
    Per-probe timeout in seconds. ``0`` lets core use its default (2s).
    """
    health_degrade_threshold: int = Field(default=3, ge=0)
    """
    Consecutive probe failures required to degrade the site. ``0`` uses the core default (3).
    """
    health_recover_threshold: int = Field(default=2, ge=0)
    """
    Consecutive probe successes required to recover the site. ``0`` uses the core default (2).
    """
    detach: Optional[bool] = Field(default=None)
    """
    Run ``swanlab-core`` as a shared detached daemon instead of binding it to this run's process lifetime.
    ``True`` (Detached Mode): runs for the same user on the same host share one ``swanlab-core`` process,
    and ``close()`` / ``abort()`` only release this client's session and channel; ``False`` (Owner Mode):
    each run spawns an isolated ``swanlab-core`` and tears it down on close.

    ``None`` (default) leaves the choice to ``run.parallel``: ``parallel="shared"`` shares one detached
    core across ranks, any other value keeps one core per process. Set ``True`` / ``False`` to override.
    """
    idle_timeout: float = Field(default=30.0, ge=0)
    """
    Idle shutdown timeout (in seconds) passed when spawning a detached ``swanlab-core``
    (applies to both explicit ``core.detach=True`` and implicit detachment via ``run.parallel="shared"``).
    ``0`` disables automatic idle shutdown.
    """
    model_config = ConfigDict(frozen=True)
