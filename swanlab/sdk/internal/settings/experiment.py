"""
@author: cunyue
@file: experiment.py
@time: 2026/3/5 20:17
@description: SwanLab 实验配置，配置这一次实验、运行的相关信息
这部分配置与业务强相关，因此基本交给context业务上下文处理，这里的职责主要是接入环境变量
"""

import json
from pathlib import Path
from typing import Annotated, Any, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator
from pydantic.config import ConfigDict
from pydantic_settings import NoDecode

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg import constraints as const
from swanlab.sdk.typings.run import ParallelType, ResumeType

from .compat import (
    experiment_color_factory,
    experiment_description_factory,
    experiment_group_factory,
    experiment_job_type_factory,
    experiment_name_factory,
    experiment_tags_factory,
    map_resume_value,
    project_name_factory,
    project_public_factory,
    run_id_factory,
    run_resume_factory,
    workspace_factory,
)


class ProjectSettings(BaseModel):
    name: Optional[const.ProjectName] = Field(
        default_factory=project_name_factory,
        validate_default=True,
    )
    """
    Project name for this SwanLab run.
    """

    workspace: Optional[const.Workspace] = Field(
        default_factory=workspace_factory,
        validate_default=True,
    )
    """
    Workspace name for this SwanLab run belongs to.
    """

    public: bool = Field(default_factory=project_public_factory)
    """
    Whether this SwanLab run is public.
    """

    model_config = ConfigDict(frozen=True)


Tags = Field(default_factory=experiment_tags_factory, max_length=30, validate_default=True)


class ExperimentSettings(BaseModel):
    name: Optional[const.ExperimentName] = Field(
        default_factory=experiment_name_factory,
        validate_default=True,
    )
    """
    Experiment name for this SwanLab run.
    """

    color: Optional[const.HexColor] = Field(
        default_factory=experiment_color_factory,
        validate_default=True,
    )
    """
    Color for this SwanLab run.
    """

    description: Optional[const.Description] = Field(
        default_factory=experiment_description_factory,
        validate_default=True,
    )
    """
    Description for this SwanLab run.
    """

    tags: Annotated[List[const.TagString], NoDecode] = Tags
    """
    Tags for this SwanLab run.
    """

    @field_validator("tags", mode="before")
    def validate_tags(cls, v: Any) -> List[str]:
        """
        自定义标签解析，同时支持JSON、逗号分隔的字符串等格式
        """
        if isinstance(v, list):
            return v
        if isinstance(v, dict):
            return [json.dumps(v)]
        # 尝试解析为JSON，并且判断是否为数组
        if isinstance(v, str):
            v = v.strip()
            if not v:
                return []
            try:
                json_value = json.loads(v)
                if isinstance(json_value, list):
                    return json_value
            except json.JSONDecodeError:
                pass
            return [item.strip() for item in v.split(",") if item.strip()]
        raise ValueError(f"tags must be a list, dict, or string, but got {type(v).__name__}")

    group: Optional[const.Group] = Field(default_factory=experiment_group_factory, validate_default=True)
    """
    Group for this SwanLab run.
    """

    job_type: Optional[const.JobType] = Field(default_factory=experiment_job_type_factory, validate_default=True)
    """
    Job type for this SwanLab run.
    """

    model_config = ConfigDict(frozen=True)


class RunSettings(BaseModel):
    id: Optional[const.RunId] = Field(default_factory=run_id_factory)
    """
    Run ID for this SwanLab run.
    """

    resume: ResumeType = Field(default_factory=run_resume_factory)
    """
    Resume policy for this SwanLab run.
    """

    @field_validator("resume", mode="before")
    def validate_resume(cls, v: Any) -> Any:
        return map_resume_value(v)

    parallel: ParallelType = Field(default="none")
    """
    Parallel execution policy for this SwanLab run.
    """

    history: Optional[Path] = Field(default=None)
    """
    History path for this SwanLab run, useful for resuming from a previous run.
    """

    config: Optional[Path] = Field(default=None)
    """
    Config file path or dict for this SwanLab run.
    """

    dir: Optional[const.RunDir] = Field(default=None)
    """
    Custom run directory name. When specified, dir_create_retries is ignored
    and the directory is created exactly once; creation failure raises an error.
    """

    dir_max_length: int = Field(default=255, ge=50, le=255)
    """
    Maximum length for the generated run directory name.
    """

    dir_create_retries: int = Field(default=10, ge=1)
    """
    Maximum number of retries for creating a unique run directory.
    If the generated directory name conflicts with an existing one, a new name will be generated after a short delay, up to this many times.
    """

    @model_validator(mode="before")
    @classmethod
    def force_resume_for_parallel(cls, data: Any) -> Any:
        if isinstance(data, dict) and data.get("parallel", "none") != "none":
            console.debug("Parallel execution detected, forcing resume policy to 'allow'")
            return {**data, "resume": "allow"}
        return data

    @model_validator(mode="after")
    def validate_dir_length(self) -> "RunSettings":
        if self.dir is not None and len(self.dir) > self.dir_max_length:
            raise ValueError(f"run.dir length ({len(self.dir)}) exceeds run.dir_max_length ({self.dir_max_length})")
        return self

    model_config = ConfigDict(frozen=True)
