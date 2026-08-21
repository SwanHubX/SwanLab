"""
实验域遗留环境变量兼容层。

0.7.x 及以前 SwanLab 使用一批与字段名不对应的环境变量（如 ``SWANLAB_PROJ_NAME`` 对应 ``project.name``），
pydantic-settings 无法映射，故通过 default_factory 直接读取 ``os.environ``。
这些读取统一下沉到本模块，并经 :func:`swanlab.sdk.internal.settings.compat.env.getenv` 剥离成对引号。
"""

from typing import Any, List, Optional, cast, get_args

from pydantic import ValidationError

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg import constraints as const
from swanlab.sdk.typings.run import ResumeType

from .env import getenv

__all__ = [
    "experiment_color_factory",
    "experiment_description_factory",
    "experiment_group_factory",
    "experiment_job_type_factory",
    "experiment_name_factory",
    "experiment_tags_factory",
    "map_resume_value",
    "project_name_factory",
    "project_public_factory",
    "run_id_factory",
    "run_resume_factory",
    "workspace_factory",
]


def project_name_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    # 这里需要额外说明：
    # 在0.7.x及以前环境变量中，即使设置了错误的 SWANLAB_PROJ_NAME ，也不会影响例如 swanlab.login 等命令
    # 这是因为曾经swanlab没有一个全局settings解析，是按需、随用随解析的——而在0.8.x及以后，swanlab对用户参数的解析放到了全局
    # 尽管这在大多数情况下不会出现行为差异，但是在一些极端情况下会出现错误
    # 其中以 SWANLAB_PROJ_NAME 最为显著，有些行为会在login前设置一个错误（例如携带双引号）的 SWANLAB_PROJ_NAME
    # 然后在login后，又手动设置回正确的 SWANLAB_PROJ_NAME —— 虽然这确实有些奇怪，但是确实导致了一些行为差异
    # 因此这里对 SWANLAB_PROJ_NAME 做了一个特殊处理：先剥离成对引号，再尝试将其解析为一个合法的项目名，解析失败则返回 None
    value = getenv("SWANLAB_PROJ_NAME", None)
    if not value:
        return None
    try:
        return const.ta_project.validate_python(value)
    except ValidationError:
        console.warning(f"Invalid SWANLAB_PROJ_NAME='{value}', ignored, using default project name instead.")
        return None


def workspace_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_WORKSPACE", None)


def project_public_factory() -> bool:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_PUBLIC", "").lower() in ["true", "yes", "1"]


def experiment_name_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_EXP_NAME", None)


def experiment_color_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_EXP_COLOR", None)


def experiment_description_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_DESCRIPTION", None)


def experiment_tags_factory() -> List[str]:
    # 向下兼容旧版本环境变量
    env_value = getenv("SWANLAB_TAGS", "")
    return [item.strip() for item in env_value.split(",") if item.strip()]


def experiment_group_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_GROUP", None)


def experiment_job_type_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_JOB_TYPE", None)


def run_id_factory() -> Optional[str]:
    # 向下兼容旧版本环境变量
    return getenv("SWANLAB_RUN_ID", None)


def run_resume_factory() -> ResumeType:
    # 向下兼容旧版本环境变量
    current_value = getenv("SWANLAB_RESUME", "").lower()
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
        if value_lower in get_args(ResumeType):
            return cast(ResumeType, value_lower)
        raise ValueError(f"Invalid resume value: {value_lower}, must be one of ['must', 'allow', 'never']")
    raise ValueError(f"Invalid resume value type: {type(value).__name__}, must be one of ['must', 'allow', 'never']")
