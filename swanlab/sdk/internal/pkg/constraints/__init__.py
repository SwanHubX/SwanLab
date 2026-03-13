"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 19:08
@description: SwanLab SDK 字段约束定义，使用 Pydantic V2 Annotated 类型声明式地定义校验规则，
作为 settings 各子模型字段约束的单一来源（Single Source of Truth）。
"""

from typing import Annotated

from pydantic import Field, TypeAdapter
from pydantic.functional_validators import AfterValidator

__all__ = [
    # Annotated types
    "ProjectName",
    "Workspace",
    "ExperimentName",
    "HexColor",
    "Description",
    "TagString",
    "Group",
    "JobType",
    "RunId",
    "MetricKey",
    "Label",
    # TypeAdapters
    "ta_project",
    "ta_workspace",
    "ta_experiment",
    "ta_description",
    "ta_hex_color",
    "ta_group",
    "ta_job_type",
    "ta_run_id",
    "ta_metric_key",
    "ta_label",
    # Validation functions
    "validate_project",
    "validate_workspace",
    "validate_experiment",
    "validate_description",
    "validate_job_type",
    "validate_group",
    "validate_run_id",
    "validate_metric_key",
]

# ---------------------------------------------------------------------------
# Annotated types — single source of truth for field constraints
# ---------------------------------------------------------------------------

ProjectName = Annotated[
    str,
    Field(min_length=1, max_length=100, pattern=r"^[0-9a-zA-Z_\-+.]+$"),
]
"""Project name: 1-100 chars, alphanumeric + ``_ - + .``."""

Workspace = Annotated[
    str,
    Field(min_length=1, max_length=25, pattern=r"^[0-9a-zA-Z\-_]+$"),
]
"""Workspace name: 1-25 chars, alphanumeric + ``- _``."""

ExperimentName = Annotated[
    str,
    Field(min_length=1, max_length=250),
]
"""Experiment name: 1-250 chars."""

HexColor = Annotated[
    str,
    Field(min_length=7, max_length=7, pattern=r"^#[0-9a-fA-F]{6}$"),
]
"""Hex color string: exactly ``#RRGGBB`` format."""

Description = Annotated[
    str,
    Field(min_length=1, max_length=1024),
]
"""Description: 1-1024 chars."""

TagString = Annotated[str, Field(max_length=200)]
"""Single tag string: up to 200 chars."""

Group = Annotated[str, Field(min_length=1, max_length=256)]
"""Experiment group: 1-256 chars."""

JobType = Annotated[str, Field(min_length=1, max_length=256)]
"""Job type: 1-256 chars."""

RunId = Annotated[
    str,
    Field(min_length=1, max_length=64, pattern=r"^[^/\\#?%:]+$"),
]
"""Run ID: 1-64 chars, must not contain ``/``, ``\\``, ``#``, ``?``, ``%``, or ``:``."""


def _no_dot_slash_edges(v: str) -> str:
    if v.startswith((".", "/")) or v.endswith((".", "/")):
        raise ValueError(f"Key '{v}' cannot start or end with '.' or '/'.")
    return v


MetricKey = Annotated[
    str,
    Field(min_length=1, max_length=255),
    AfterValidator(_no_dot_slash_edges),
]
"""Metric / log key: 1-255 chars, must not start or end with ``'.'`` or ``'/'``."""

Label = Annotated[str, Field(max_length=255)]
"""Display label for metrics and charts: up to 255 chars."""

# ---------------------------------------------------------------------------
# TypeAdapters (module-level singletons — build once, reuse everywhere)
# ---------------------------------------------------------------------------

ta_project = TypeAdapter(ProjectName)
ta_workspace = TypeAdapter(Workspace)
ta_experiment = TypeAdapter(ExperimentName)
ta_description = TypeAdapter(Description)
ta_hex_color = TypeAdapter(HexColor)
ta_group = TypeAdapter(Group)
ta_job_type = TypeAdapter(JobType)
ta_run_id = TypeAdapter(RunId)
ta_metric_key = TypeAdapter(MetricKey)
ta_label = TypeAdapter(Label)

# ---------------------------------------------------------------------------
# Validation helpers — thin wrappers around TypeAdapter.validate_python()
# ---------------------------------------------------------------------------


def validate_project(name: str) -> str:
    return ta_project.validate_python(name.strip() if isinstance(name, str) else name)


def validate_workspace(workspace: str) -> str:
    return ta_workspace.validate_python(workspace.strip() if isinstance(workspace, str) else workspace)


def validate_experiment(name: str) -> str:
    return ta_experiment.validate_python(name.strip() if isinstance(name, str) else name)


def validate_description(desc: str) -> str:
    return ta_description.validate_python(desc)


def validate_group(group: str) -> str:
    return ta_group.validate_python(group)


def validate_job_type(job_type: str) -> str:
    return ta_job_type.validate_python(job_type)


def validate_run_id(run_id: str) -> str:
    return ta_run_id.validate_python(run_id.strip() if isinstance(run_id, str) else run_id)


def validate_metric_key(key: str) -> str:
    return ta_metric_key.validate_python(key.strip() if isinstance(key, str) else key)
