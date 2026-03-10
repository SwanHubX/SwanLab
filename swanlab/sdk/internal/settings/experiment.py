"""
@author: cunyue
@file: experiment.py
@time: 2026/3/5 20:17
@description: SwanLab 实验配置，配置这一次实验、运行的相关信息
这部分配置与业务强相关，因此基本交给context业务上下文处理，这里的职责主要是接入环境变量
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, List, Literal, Optional, cast

if sys.version_info >= (3, 9):
    from typing import Annotated
else:
    from typing_extensions import Annotated

from pydantic import BaseModel, Field, field_validator
from pydantic_settings import NoDecode

from swanlab.sdk.typings.run import ResumeType


def project_name_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_PROJ_NAME", None)


def workspace_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_WORKSPACE", None)


def project_public_factory() -> bool:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_PUBLIC", "").lower() in ["true", "yes", "1"]


class ProjectSettings(BaseModel):
    name: Optional[str] = Field(
        default_factory=project_name_factory,
        max_length=100,
        pattern=r"^[0-9a-zA-Z_\-+.]+$",
        validate_default=True,
    )
    """
    Project name for this SwanLab run.
    """

    workspace: Optional[str] = Field(default_factory=workspace_factory)
    """
    Workspace name for this SwanLab run belongs to.
    """

    public: bool = Field(default_factory=project_public_factory)
    """
    Whether this SwanLab run is public.
    """


def experiment_name_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_EXP_NAME", None)


def experiment_color_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_EXP_COLOR", None)


def experiment_description_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_DESCRIPTION", None)


def experiment_tags_factory() -> List[str]:
    # 向下兼容旧版本环境变量
    env_value = os.environ.get("SWANLAB_TAGS", "")
    return [item.strip() for item in env_value.split(",") if item.strip()]


def experiment_group_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_GROUP", None)


def experiment_job_type_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_JOB_TYPE", None)


ValidTagString = Annotated[str, Field(max_length=200)]

ValidTags = Field(default_factory=experiment_tags_factory, max_length=50, validate_default=True)


class ExperimentSettings(BaseModel):
    name: Optional[str] = Field(
        default_factory=experiment_name_factory,
        min_length=1,
        max_length=250,
        validate_default=True,
    )
    """
    Experiment name for this SwanLab run.
    """

    color: Optional[str] = Field(
        default_factory=experiment_color_factory,
        max_length=7,
        min_length=7,
        pattern=r"^#[0-9a-fA-F]{6}$",
        validate_default=True,
    )
    """
    Color for this SwanLab run.
    """

    description: Optional[str] = Field(
        default_factory=experiment_description_factory,
        min_length=1,
        max_length=1024,
        validate_default=True,
    )
    """
    Description for this SwanLab run.
    """
    tags: Annotated[List[ValidTagString], NoDecode] = ValidTags
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

    group: Optional[str] = Field(
        default_factory=experiment_group_factory, min_length=1, max_length=256, validate_default=True
    )
    """
    Group for this SwanLab run.
    """

    job_type: Optional[str] = Field(
        default_factory=experiment_job_type_factory, min_length=1, max_length=256, validate_default=True
    )
    """
    Job type for this SwanLab run.
    """


def run_id_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_RUN_ID", "")


def run_resume_factory() -> ResumeType:
    # 向下兼容旧版本环境变量
    current_value = os.environ.get("SWANLAB_RESUME", "").lower()
    if current_value:
        try:
            return map_resume_value(current_value)
        except ValueError:
            pass
    return "never"


def map_resume_value(value: Any) -> ResumeType:
    """
    将各种形式的 resume 值映射为允许的 ResumeType
    """
    if isinstance(value, bool):
        return "allow" if value else "never"
    if isinstance(value, str):
        value_lower = value.lower()
        if value_lower == "true":
            return "allow"
        if value_lower == "false":
            return "never"
        if value_lower == "yes":
            return "allow"
        if value_lower == "no":
            return "never"
        if value_lower == "1":
            return "allow"
        if value_lower == "0":
            return "never"
        if value_lower in ["must", "allow", "never"]:
            return cast(ResumeType, value_lower)
        raise ValueError(f"Invalid resume value: {value_lower}, must be one of ['must', 'allow', 'never']")
    raise ValueError(f"Invalid resume value type: {type(value).__name__}, must be one of ['must', 'allow', 'never']")


class RunSettings(BaseModel):
    id: Optional[str] = Field(default_factory=run_id_factory)
    """
    Run ID for this SwanLab run.
    """

    resume: Literal["must", "allow", "never"] = Field(default_factory=run_resume_factory)
    """
    Resume policy for this SwanLab run.
    """

    @field_validator("resume", mode="before")
    def validate_resume(cls, v: Any) -> Any:
        return map_resume_value(v)

    config: Optional[Path] = Field(default=None)
