"""
@author: cunyue
@file: test_constraints.py
@time: 2026/3/13
@description: 测试 swanlab.sdk.internal.pkg.constraints 字段约束模块
"""

import pytest
from pydantic import TypeAdapter, ValidationError

from swanlab.sdk.internal.pkg.constraints import (
    Description,
    ExperimentName,
    Group,
    HexColor,
    JobType,
    MetricKey,
    ProjectName,
    RunId,
    TagString,
    Workspace,
    validate_description,
    validate_experiment,
    validate_group,
    validate_job_type,
    validate_metric_key,
    validate_project,
    validate_run_id,
    validate_workspace,
)

# ---------------------------------------------------------------------------
# ProjectName
# ---------------------------------------------------------------------------


class TestProjectName:
    adapter = TypeAdapter(ProjectName)

    @pytest.mark.parametrize(
        "value",
        ["my-project", "Project_1", "foo.bar+baz", "A" * 100, "abc123"],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "A" * 101,  # too long
            "invalid name",  # space
            "name@invalid",  # @
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_project_strips_whitespace(self):
        assert validate_project("  myproject  ") == "myproject"

    def test_validate_project_valid(self):
        assert validate_project("valid-project") == "valid-project"

    def test_validate_project_invalid(self):
        with pytest.raises(ValidationError):
            validate_project("invalid name!")


# ---------------------------------------------------------------------------
# Workspace
# ---------------------------------------------------------------------------


class TestWorkspace:
    adapter = TypeAdapter(Workspace)

    @pytest.mark.parametrize(
        "value",
        ["my-workspace", "workspace_1", "Workspace", "a" * 25],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "a" * 26,  # too long
            "invalid name",  # space
            "invalid+name",  # +
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_workspace_strips_whitespace(self):
        assert validate_workspace("  myspace  ") == "myspace"


# ---------------------------------------------------------------------------
# ExperimentName
# ---------------------------------------------------------------------------


class TestExperimentName:
    adapter = TypeAdapter(ExperimentName)

    @pytest.mark.parametrize(
        "value",
        ["experiment-1", "My Experiment", "实验名称", "x" * 250],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "x" * 251,  # too long
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_experiment_strips_whitespace(self):
        assert validate_experiment("  my experiment  ") == "my experiment"


# ---------------------------------------------------------------------------
# HexColor
# ---------------------------------------------------------------------------


class TestHexColor:
    adapter = TypeAdapter(HexColor)

    @pytest.mark.parametrize(
        "value",
        ["#528d59", "#FFFFFF", "#000000", "#aAbBcC", "#1a2b3c"],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "#FFF",  # too short
            "#GGGGGG",  # invalid hex
            "528d59",  # missing #
            "#528d590",  # too long
            "",
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)


# ---------------------------------------------------------------------------
# Description
# ---------------------------------------------------------------------------


class TestDescription:
    adapter = TypeAdapter(Description)

    def test_valid(self):
        assert self.adapter.validate_python("a description") == "a description"
        assert self.adapter.validate_python("x" * 1024) == "x" * 1024

    def test_too_long(self):
        with pytest.raises(ValidationError):
            self.adapter.validate_python("x" * 1025)

    def test_empty_invalid(self):
        with pytest.raises(ValidationError):
            self.adapter.validate_python("")

    def test_validate_description_valid(self):
        assert validate_description("hello") == "hello"

    def test_validate_description_too_long(self):
        with pytest.raises(ValidationError):
            validate_description("x" * 1025)


# ---------------------------------------------------------------------------
# TagString / Tags
# ---------------------------------------------------------------------------


class TestTagString:
    adapter = TypeAdapter(TagString)

    def test_valid(self):
        assert self.adapter.validate_python("tag") == "tag"
        assert self.adapter.validate_python("t" * 200) == "t" * 200

    def test_too_long(self):
        with pytest.raises(ValidationError):
            self.adapter.validate_python("t" * 201)


# ---------------------------------------------------------------------------
# Group
# ---------------------------------------------------------------------------


class TestGroup:
    adapter = TypeAdapter(Group)

    def test_valid(self):
        assert self.adapter.validate_python("my-group") == "my-group"
        assert self.adapter.validate_python("g" * 256) == "g" * 256

    @pytest.mark.parametrize("value", ["", "g" * 257])
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_group(self):
        assert validate_group("experiment-group") == "experiment-group"


# ---------------------------------------------------------------------------
# JobType
# ---------------------------------------------------------------------------


class TestJobType:
    adapter = TypeAdapter(JobType)

    def test_valid(self):
        assert self.adapter.validate_python("training") == "training"
        assert self.adapter.validate_python("j" * 256) == "j" * 256

    @pytest.mark.parametrize("value", ["", "j" * 257])
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_job_type(self):
        assert validate_job_type("eval") == "eval"


# ---------------------------------------------------------------------------
# RunId
# ---------------------------------------------------------------------------


class TestRunId:
    adapter = TypeAdapter(RunId)

    @pytest.mark.parametrize(
        "value",
        ["abc123", "run_001", "a" * 64, "MyRun-1", "run with spaces", "实验001"],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "a" * 65,  # too long
            "run/id",  # contains /
            "run\\id",  # contains \
            "run#id",  # contains #
            "run?id",  # contains ?
            "run%id",  # contains %
            "run:id",  # contains :
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_run_id_strips_whitespace(self):
        assert validate_run_id("  abc123  ") == "abc123"

    def test_validate_run_id_valid(self):
        assert validate_run_id("run_001") == "run_001"

    def test_validate_run_id_uppercase_valid(self):
        assert validate_run_id("RunID") == "RunID"

    def test_validate_run_id_forbidden_chars(self):
        for char in ["/", "\\", "#", "?", "%", ":"]:
            with pytest.raises(ValidationError):
                validate_run_id(f"run{char}id")


# ---------------------------------------------------------------------------
# MetricKey
# ---------------------------------------------------------------------------


class TestMetricKey:
    adapter = TypeAdapter(MetricKey)

    @pytest.mark.parametrize(
        "value",
        ["loss", "train/loss", "metrics.accuracy", "a/b.c", "x" * 255],
    )
    def test_valid(self, value):
        assert self.adapter.validate_python(value) == value

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "x" * 256,  # too long
            ".starts_with_dot",
            "/starts_with_slash",
            "ends_with_dot.",
            "ends_with_slash/",
        ],
    )
    def test_invalid(self, value):
        with pytest.raises(ValidationError):
            self.adapter.validate_python(value)

    def test_validate_metric_key_strips_whitespace(self):
        assert validate_metric_key("  loss  ") == "loss"

    def test_validate_metric_key_valid(self):
        assert validate_metric_key("train/loss") == "train/loss"

    def test_validate_metric_key_edge_dots(self):
        with pytest.raises(ValidationError):
            validate_metric_key(".hidden")
        with pytest.raises(ValidationError):
            validate_metric_key("hidden.")

    def test_validate_metric_key_edge_slashes(self):
        with pytest.raises(ValidationError):
            validate_metric_key("/root")
        with pytest.raises(ValidationError):
            validate_metric_key("root/")


# ---------------------------------------------------------------------------
# Integration: Settings models use constraints types
# ---------------------------------------------------------------------------


class TestSettingsUsesConstraints:
    """
    验证 Settings 子模型的字段实际由 constraints 中的 Annotated 类型驱动。
    """

    def test_project_name_constraint_respected(self):
        from swanlab.sdk.internal.settings.experiment import ProjectSettings

        with pytest.raises(ValidationError):
            ProjectSettings(name="invalid name!")  # space not allowed

    def test_project_name_valid(self):
        from swanlab.sdk.internal.settings.experiment import ProjectSettings

        s = ProjectSettings(name="valid-project")
        assert s.name == "valid-project"

    def test_experiment_name_too_long(self):
        from swanlab.sdk.internal.settings.experiment import ExperimentSettings

        with pytest.raises(ValidationError):
            ExperimentSettings(name="x" * 251)

    def test_run_id_constraint_respected(self):
        from swanlab.sdk.internal.settings.experiment import RunSettings

        with pytest.raises(ValidationError):
            RunSettings(id="run/id")  # / not allowed

    def test_run_id_valid(self):
        from swanlab.sdk.internal.settings.experiment import RunSettings

        s = RunSettings(id="My-Run_001")
        assert s.id == "My-Run_001"

    def test_hex_color_constraint_respected(self):
        from swanlab.sdk.internal.settings.experiment import ExperimentSettings

        with pytest.raises(ValidationError):
            ExperimentSettings(color="#ZZZ")

    def test_hex_color_valid(self):
        from swanlab.sdk.internal.settings.experiment import ExperimentSettings

        s = ExperimentSettings(color="#528d59")
        assert s.color == "#528d59"
