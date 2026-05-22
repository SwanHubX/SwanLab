"""
@author: cunyue
@file: __init__.py
@time: 2026/5/11 19:19
@description: 探针上下文
"""

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Optional

from swanlab.proto.swanlab.settings.probe.v1.probe_pb2 import ProbeSettings
from swanlab.sdk.internal.pkg import console

__all__ = ["ProbeConfig", "ProbeContext"]


@dataclass(frozen=True)
class ProbeConfig:
    run_id: Optional[str]
    run_dir: Optional[Path]
    global_system_step: int
    hardware: bool
    runtime: bool
    requirements: bool
    conda: bool
    git: bool
    swanlab: bool
    monitor: bool
    monitor_interval: int
    monitor_disk_dir: Path


class ProbeContext:
    def __init__(self, *, config: ProbeConfig):
        self.config = config

    @classmethod
    def from_proto(cls, proto: ProbeSettings) -> "ProbeContext":
        run_id = proto.run_id if proto.HasField("run_id") else None
        run_dir = Path(proto.run_dir) if proto.HasField("run_dir") else None
        if proto.HasField("global_system_step"):
            global_system_step = proto.global_system_step
        else:
            console.debug("ProbeSettings missing global_system_step, defaulting to 0")
            global_system_step = 0

        config = ProbeConfig(
            run_id=run_id,
            run_dir=run_dir,
            global_system_step=global_system_step,
            hardware=proto.hardware,
            runtime=proto.runtime,
            requirements=proto.requirements,
            conda=proto.conda,
            git=proto.git,
            swanlab=proto.swanlab,
            monitor=proto.monitor,
            monitor_interval=proto.monitor_interval,
            monitor_disk_dir=Path(proto.monitor_disk_dir),
        )
        return cls(config=config)

    @cached_property
    def files_dir(self) -> Path:
        assert self.config.run_dir is not None, "run_dir must be set in ProbeConfig to access files_dir"
        return self.config.run_dir / "files"

    @cached_property
    def metadata_file(self) -> Path:
        return self.files_dir / "swanlab-metadata.json"

    @cached_property
    def requirements_file(self) -> Path:
        return self.files_dir / "requirements.txt"

    @cached_property
    def conda_file(self) -> Path:
        return self.files_dir / "conda.yaml"
