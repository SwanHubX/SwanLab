"""
@author: cunyue
@file: __init__.py
@time: 2026/5/11 18:19
@description: core 运行上下文
"""

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings

__all__ = ["CoreContext", "CoreConfig"]


@dataclass(frozen=True)
class CoreConfig:
    run_id: str
    run_dir: Path
    section_rule: int
    record_batch: int
    record_interval: float
    save_size: int
    save_split: int
    save_part: int
    save_batch: int


class CoreContext:
    def __init__(self, *, config: CoreConfig):
        self.config = config

    @classmethod
    def from_proto(cls, proto: CoreSettings) -> "CoreContext":
        config = CoreConfig(
            run_id=proto.run_id,
            run_dir=Path(proto.run_dir),
            section_rule=proto.section_rule,
            record_batch=proto.record_batch,
            record_interval=proto.record_interval,
            save_size=proto.save_size,
            save_split=proto.save_split,
            save_part=proto.save_part,
            save_batch=proto.save_batch,
        )
        return cls(config=config)

    @cached_property
    def media_dir(self) -> Path:
        return self.config.run_dir / "media"

    @cached_property
    def debug_dir(self) -> Path:
        return self.config.run_dir / "debug"

    @cached_property
    def files_dir(self) -> Path:
        return self.config.run_dir / "files"

    @cached_property
    def metadata_file(self) -> Path:
        return self.files_dir / "swanlab-metadata.json"

    @cached_property
    def config_file(self) -> Path:
        return self.files_dir / "config.yaml"

    @cached_property
    def requirements_file(self) -> Path:
        return self.files_dir / "requirements.txt"

    @cached_property
    def conda_file(self) -> Path:
        return self.files_dir / "conda.yaml"

    @cached_property
    def run_file(self) -> Path:
        assert self.config.run_id, "Run ID is not set."
        return self.config.run_dir / f"run-{self.config.run_id}.swanlab"
