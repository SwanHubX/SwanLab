"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 19:08
@description: SwanLab SDK 字段约束定义，使用 Pydantic V2 Annotated 类型声明式地定义校验规则，
作为 settings 各子模型字段约束的单一来源（Single Source of Truth）。
"""

import re
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
    # MetricKey sanitization regex
    "METRIC_KEY_INVALID_RE",
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
    "ta_chart_name",
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


# 合法字符集：word chars（\w）、点、斜杠、连字符
# 同时作为 MetricKey pattern 约束和 sanitize 用的反向正则的来源
_METRIC_KEY_VALID_CHARS = r"\w./-"

METRIC_KEY_INVALID_RE = re.compile(rf"[^{_METRIC_KEY_VALID_CHARS}]")
"""预编译的非法字符正则，供 sanitize 场景（如 validate_key）使用。"""

MetricKey = Annotated[
    str,
    Field(min_length=1, max_length=255, pattern=rf"^[{_METRIC_KEY_VALID_CHARS}]+$"),
    AfterValidator(_no_dot_slash_edges),
]
"""Metric / log key: 1-255 chars, must not start or end with ``'.'`` or ``'/'``."""

Label = Annotated[str, Field(max_length=255)]
"""Display label for metrics and runs: up to 255 chars."""

ChartName = Annotated[str, Field(min_length=1, max_length=255)]
"""Display name for charts: 1-255 chars."""

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
ta_chart_name = TypeAdapter(ChartName)
