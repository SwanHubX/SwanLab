"""
@author: cunyue
@file: __init__.py
@time: 2026/5/11 18:19
@description: core 运行上下文
"""

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Literal, Optional

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
    def __init__(self, *, config: CoreConfig, mode: Literal["core", "sync"]):
        self.config = config
        self._mode: Literal["core", "sync"] = mode
        # 云端信息
        self._username: Optional[str] = None
        self._project: Optional[str] = None
        self._project_id: Optional[str] = None
        self._experiment_id: Optional[str] = None

    @classmethod
    def from_proto(cls, proto: CoreSettings, mode: Literal["core", "sync"] = "core") -> "CoreContext":
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
        return cls(config=config, mode=mode)

    def set_online_params(self, username: str, project: str, project_id: str, experiment_id: str):
        """
        设置云端信息
        :param username: 用户名
        :param project: 项目名
        :param project_id: 项目id
        :param experiment_id: 实验id
        """
        self._username = username
        self._project = project
        self._project_id = project_id
        self._experiment_id = experiment_id

    @cached_property
    def username(self) -> str:
        assert self._username is not None, (
            "ctx.username is not initialized. "
            "This usually means the current run is not in online mode, "
            "or an internal initialization bug occurred."
        )
        return self._username

    @cached_property
    def project(self) -> str:
        assert self._project is not None, (
            "ctx.project is not initialized. "
            "This usually means the current run is not in online mode, "
            "or an internal initialization bug occurred."
        )
        return self._project

    @cached_property
    def project_id(self) -> str:
        assert self._project_id is not None, (
            "ctx.project_id is not initialized. "
            "This usually means the current run is not in online mode, "
            "or an internal initialization bug occurred."
        )
        return self._project_id

    @cached_property
    def experiment_id(self) -> str:
        assert self._experiment_id is not None, (
            "ctx.experiment_id is not initialized. "
            "This usually means the current run is not in online mode, "
            "or an internal initialization bug occurred."
        )
        return self._experiment_id

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
        if self._mode == "core":
            assert self.config.run_id, "Run ID is not set."
            # core 模式下，根据run id创建run file
            return self.config.run_dir / f"run-{self.config.run_id}.swanlab"
        elif self._mode == "sync":
            # sync 模式下，查询目录下 run-*.swanlab 文件作为run file
            files = sorted(self.config.run_dir.glob("run-*.swanlab"))
            if len(files) == 0:
                raise FileNotFoundError(f"No run-*.swanlab file found in {self.config.run_dir}")
            if len(files) > 1:
                raise RuntimeError(f"Multiple run-*.swanlab files found in {self.config.run_dir}")
            return files[0]
        raise ValueError(f"Unsupported CoreContext mode: {self._mode}")
