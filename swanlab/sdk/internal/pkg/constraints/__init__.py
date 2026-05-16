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
    "RunDir",
    "MetricKey",
    "MetricName",
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
    "ta_metric_name",
    "ta_chart_name",
    # value
    "METRIC_KEY_MAX_LENGTH",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
METRIC_KEY_MAX_LENGTH = 512

OS_SAFE_PATTERN = r'^[^<>:"/\\|?*#%\x00-\x1f\x7f]+$'

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
    Field(min_length=1, max_length=512),
]
"""Experiment name: 1-512 chars."""

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

TagString = Annotated[str, Field(max_length=20)]
"""Single tag string: up to 20 chars."""

Group = Annotated[str, Field(min_length=1, max_length=512)]
"""Experiment group: 1-512 chars."""

JobType = Annotated[str, Field(min_length=1, max_length=512)]
"""Job type: 1-512 chars."""

RunId = Annotated[
    str,
    Field(min_length=1, max_length=512, pattern=OS_SAFE_PATTERN),
]
"""Run ID: 1-512 chars, must not contain OS-reserved characters or control characters."""

RunDir = Annotated[
    str,
    Field(min_length=1, max_length=255, pattern=OS_SAFE_PATTERN),
]
"""Custom run directory name: 1-255 chars, must not contain OS-reserved characters or control characters."""


def _no_dot_slash_edges(v: str) -> str:
    if v.startswith((".", "/")) or v.endswith((".", "/")):
        raise ValueError(f"Key '{v}' cannot start or end with '.' or '/'.")
    return v


MetricKey = Annotated[
    str,
    Field(min_length=1, max_length=METRIC_KEY_MAX_LENGTH, pattern=r"^[^\x00-\x1f\x7f]+$"),
    AfterValidator(_no_dot_slash_edges),
]
"""Metric / log key: 1-512 chars, no control characters, must not start or end with ``'.'`` or ``'/'``."""

MetricName = Annotated[str, Field(max_length=512)]
"""Display label for metrics and runs: up to 512 chars."""

ChartName = Annotated[str, Field(min_length=1, max_length=512)]
"""Display name for charts: 1-512 chars."""

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
ta_metric_name = TypeAdapter(MetricName)
ta_chart_name = TypeAdapter(ChartName)
