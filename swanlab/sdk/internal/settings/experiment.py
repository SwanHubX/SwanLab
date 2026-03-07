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
from typing import Any, List, Literal, cast

if sys.version_info >= (3, 9):
    from typing import Annotated
else:
    from typing_extensions import Annotated

from pydantic import BaseModel, Field, field_validator
from pydantic_settings import NoDecode

from swanlab.sdk.typing.run import ResumeType


def project_name_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_PROJ_NAME", "")


def workspace_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_WORKSPACE", "")


class ProjectSettings(BaseModel):
    name: str = Field(default_factory=project_name_factory)
    """
    Project name for this SwanLab run.
    """

    workspace: str = Field(default_factory=workspace_factory)
    """
    Workspace name for this SwanLab run belongs to.
    """


def experiment_name_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_EXP_NAME", "")


def experiment_description_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_DESCRIPTION", "")


def experiment_tags_factory() -> List[str]:
    # 向下兼容旧版本环境变量
    env_value = os.environ.get("SWANLAB_TAGS", "")
    return [item.strip() for item in env_value.split(",") if item.strip()]


def experiment_group_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_GROUP", "")


def experiment_job_type_factory() -> str:
    # 向下兼容旧版本环境变量
    return os.environ.get("SWANLAB_JOB_TYPE", "")


class ExperimentSettings(BaseModel):
    name: str = Field(default_factory=experiment_name_factory)
    """
    Experiment name for this SwanLab run.
    """

    description: str = Field(default_factory=experiment_description_factory)
    """
    Description for this SwanLab run.
    """

    tags: Annotated[List[str], NoDecode] = Field(default_factory=experiment_tags_factory)  # noqa: cannot specify `Annotated` and value `Field`s together for 'tags'
    """
    Tags for this SwanLab run.
    """

    @field_validator("tags", mode="before")
    def validate_tags(cls, v: Any) -> list[str]:
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

    group: str = Field(default_factory=experiment_group_factory)
    """
    Group for this SwanLab run.
    """

    job_type: str = Field(default_factory=experiment_job_type_factory)
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
    id: str = Field(default_factory=run_id_factory)
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
