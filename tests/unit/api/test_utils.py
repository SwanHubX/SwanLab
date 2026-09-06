"""
@author: caddiesnew
@time: 2026/4/27
@description: swanlab/api 校验函数单测
"""

from typing import cast

import pytest

from swanlab.api.helper import (
    RELATIVE_TIME_AXIS,
    STEP_AXIS,
    TIME_AXIS,
    builtin_x_axis,
    parse_timestamp_ms,
    validate_column_params,
    validate_filter,
    validate_group,
    validate_metric_log_level,
    validate_metric_type,
    validate_project_name,
    validate_sort,
    validate_x_axis,
)
from swanlab.api.self_hosted import SelfHosted
from swanlab.api.typings.common import PaginatedQuery, RangeQuery
from swanlab.api.typings.selfhosted import ApiSelfHostedInfoType


# ---------------------------------------------------------------------------
# validate_project_name
# ---------------------------------------------------------------------------
class TestValidateProjectName:
    def test_valid(self):
        validate_project_name("my-project_1.0+beta")

    @pytest.mark.parametrize("name", ["", "x" * 101])
    def test_length_invalid(self, name):
        with pytest.raises(ValueError, match="1 and 100"):
            validate_project_name(name)

    @pytest.mark.parametrize("name", ["hello world", "中文项目", "a/b", "a@b"])
    def test_invalid_chars(self, name):
        with pytest.raises(ValueError, match="0-9"):
            validate_project_name(name)


# ---------------------------------------------------------------------------
# validate_column_params
# ---------------------------------------------------------------------------
class TestValidateColumnParams:
    def test_valid_type_and_class(self):
        validate_column_params(column_type="FLOAT", column_class="CUSTOM")

    def test_invalid_type(self):
        with pytest.raises(ValueError, match="Invalid column_type"):
            validate_column_params(column_type="INVALID")

    def test_invalid_class(self):
        with pytest.raises(ValueError, match="Invalid column_class"):
            validate_column_params(column_class="INVALID")


# ---------------------------------------------------------------------------
# validate_metric_type / validate_metric_log_level
# ---------------------------------------------------------------------------
class TestValidateMetricType:
    def test_valid_scalar(self):
        validate_metric_type("SCALAR", key="loss")

    def test_log_no_key_ok(self):
        validate_metric_type("LOG")

    def test_scalar_without_key_raises(self):
        with pytest.raises(ValueError, match="key is required"):
            validate_metric_type("SCALAR", key="")

    def test_invalid_type(self):
        with pytest.raises(ValueError, match="Invalid metric_type"):
            validate_metric_type("INVALID", key="x")


class TestValidateMetricLogLevel:
    def test_valid(self):
        validate_metric_log_level("INFO")

    def test_invalid(self):
        with pytest.raises(ValueError, match="Invalid metric log level"):
            validate_metric_log_level("VERBOSE")


# ---------------------------------------------------------------------------
# validate_filter / validate_group / validate_sort
# ---------------------------------------------------------------------------
class TestValidateFilter:
    def test_valid(self):
        validate_filter({"key": "name", "type": "STABLE", "op": "EQ", "value": ["test"]})

    def test_missing_fields(self):
        with pytest.raises(ValueError, match="Missing required"):
            validate_filter({"key": "name"})

    def test_invalid_type(self):
        with pytest.raises(ValueError, match="Invalid type"):
            validate_filter({"key": "name", "type": "INVALID", "op": "EQ", "value": ["x"]})

    def test_invalid_op(self):
        with pytest.raises(ValueError, match="Invalid filter op"):
            validate_filter({"key": "name", "type": "STABLE", "op": "LIKE", "value": ["x"]})

    def test_value_not_list(self):
        with pytest.raises(ValueError, match="must be a list"):
            validate_filter({"key": "name", "type": "STABLE", "op": "EQ", "value": "not_list"})

    def test_invalid_stable_key(self):
        with pytest.raises(ValueError, match="Invalid STABLE key"):
            validate_filter({"key": "invalid_key", "type": "STABLE", "op": "EQ", "value": ["x"]})


class TestValidateGroup:
    def test_valid(self):
        validate_group({"key": "cluster", "type": "STABLE"})

    def test_missing_fields(self):
        with pytest.raises(ValueError, match="Missing required"):
            validate_group({"key": "name"})


class TestValidateSort:
    def test_valid(self):
        validate_sort({"key": "name", "type": "STABLE", "order": "ASC"})

    def test_invalid_order(self):
        with pytest.raises(ValueError, match="Invalid sort order"):
            validate_sort({"key": "name", "type": "STABLE", "order": "RANDOM"})


# ---------------------------------------------------------------------------
# PaginatedQuery
# ---------------------------------------------------------------------------
class TestPaginatedQuery:
    def test_valid_defaults(self):
        q = PaginatedQuery()
        assert q.page == 1 and q.size == 100

    def test_page_less_than_1(self):
        with pytest.raises(ValueError, match="page must be >= 1"):
            PaginatedQuery(page=0)

    def test_invalid_size(self):
        with pytest.raises(ValueError, match="size must be one of"):
            PaginatedQuery(size=42)

    def test_to_params_filters_none(self):
        q = PaginatedQuery()
        params = q.to_params(search=None, sort=None)
        assert "search" not in params
        assert "sort" not in params

    def test_to_params_includes_extras(self):
        q = PaginatedQuery()
        params = q.to_params(detail=True, extra_key="val")
        assert params["detail"] is True
        assert params["extra_key"] == "val"


# ---------------------------------------------------------------------------
# SelfHosted validation
# ---------------------------------------------------------------------------
class TestSelfHostedValidation:
    def _make_info(self, **overrides) -> ApiSelfHostedInfoType:
        base: dict = {
            "enabled": True,
            "expired": False,
            "root": True,
            "plan": "free",
            "seats": 10,
        }
        base.update(overrides)
        return cast(ApiSelfHostedInfoType, base)

    def test_validate_expire_ok(self):
        SelfHosted.validate_expire(self._make_info(expired=False))

    def test_validate_expire_raises(self):
        with pytest.raises(ValueError, match="expired"):
            SelfHosted.validate_expire(self._make_info(expired=True))

    def test_validate_root_ok(self):
        SelfHosted.validate_root(self._make_info(expired=False, root=True))

    def test_validate_root_not_root(self):
        with pytest.raises(ValueError, match="root"):
            SelfHosted.validate_root(self._make_info(expired=False, root=False))

    def test_validate_root_expired(self):
        with pytest.raises(ValueError, match="expired"):
            SelfHosted.validate_root(self._make_info(expired=True, root=True))


# ---------------------------------------------------------------------------
# parse_timestamp_ms
# ---------------------------------------------------------------------------
class TestParseTimestampMs:
    def test_int_timestamp_ms(self):
        assert parse_timestamp_ms(1715769600123) == 1715769600123

    def test_int_timestamp_seconds_to_ms(self):
        assert parse_timestamp_ms(1715769600) == 1715769600_000

    def test_numeric_string_ms(self):
        assert parse_timestamp_ms("1715769600123") == 1715769600123

    def test_numeric_string_seconds_to_ms(self):
        assert parse_timestamp_ms("1715769600") == 1715769600_000

    def test_invalid_string_raises(self):
        with pytest.raises(ValueError, match="Invalid timestamp"):
            parse_timestamp_ms("not-a-timestamp")

    def test_unsupported_type_raises(self):
        with pytest.raises(ValueError, match="Expected str or int"):
            parse_timestamp_ms(None)  # type: ignore[arg-type]

    def test_float_raises(self):
        with pytest.raises(ValueError, match="Expected str or int"):
            parse_timestamp_ms(1715769600.0)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# validate_x_axis
# ---------------------------------------------------------------------------
class TestValidateXAxis:
    @pytest.mark.parametrize("axis", ["step", "time", "relative_time"])
    def test_valid_builtins(self, axis: str):
        assert validate_x_axis(axis) == axis

    @pytest.mark.parametrize("axis", ["epoch", "iter", "custom_axis", "123", "lr/step"])
    def test_valid_custom_key(self, axis: str):
        assert validate_x_axis(axis) == axis

    @pytest.mark.parametrize("axis", ["", "   ", "\t\n"])
    def test_reject_empty_or_whitespace(self, axis: str):
        with pytest.raises(ValueError, match="non-empty string"):
            validate_x_axis(axis)

    @pytest.mark.parametrize("axis", [None, 123, ["epoch"]])
    def test_reject_non_string(self, axis):
        with pytest.raises(ValueError, match="non-empty string"):
            validate_x_axis(axis)  # type: ignore[arg-type]

    @pytest.mark.parametrize("metric_type", ["MEDIA", "LOG"])
    def test_non_scalar_accepts_step(self, metric_type: str):
        assert validate_x_axis("step", metric_type=metric_type) == "step"

    @pytest.mark.parametrize("metric_type", ["MEDIA", "LOG"])
    @pytest.mark.parametrize("axis", ["time", "relative_time"])
    def test_non_scalar_rejects_builtin_time_axes(self, metric_type: str, axis: str):
        with pytest.raises(ValueError, match="only supported for SCALAR metrics"):
            validate_x_axis(axis, metric_type=metric_type)

    @pytest.mark.parametrize("metric_type", ["MEDIA", "LOG"])
    def test_non_scalar_rejects_custom_key(self, metric_type: str):
        with pytest.raises(ValueError, match="custom x_axis key 'epoch' is only supported for SCALAR metrics"):
            validate_x_axis("epoch", metric_type=metric_type)


# ---------------------------------------------------------------------------
# builtin_x_axis & axis_request_params
# ---------------------------------------------------------------------------
class TestBuiltinXAxis:
    def test_builtins(self):
        assert builtin_x_axis("step") == STEP_AXIS
        assert builtin_x_axis("time") == TIME_AXIS
        assert builtin_x_axis("relative_time") == RELATIVE_TIME_AXIS

    @pytest.mark.parametrize("axis", ["epoch", "auto", "timestamp", "unknown", ""])
    def test_non_builtins_return_none(self, axis: str):
        assert builtin_x_axis(axis) is None


# ---------------------------------------------------------------------------
# RangeQuery custom & bounds validation
# ---------------------------------------------------------------------------
class TestRangeQueryBounds:
    def test_step_integer_bounds_accepted(self):
        rq = RangeQuery(type="step", start=0, end=100)
        assert rq.start == 0.0
        assert rq.end == 100.0

        rq_float_int = RangeQuery(type="step", start=5.0, end=10.0)
        assert rq_float_int.start == 5.0
        assert rq_float_int.end == 10.0

    @pytest.mark.parametrize("invalid_val", [1.5, -1, -0.1])
    def test_step_rejects_fraction_or_negative(self, invalid_val):
        with pytest.raises(ValueError, match="must be a non-negative integer for type 'step'"):
            RangeQuery(type="step", start=invalid_val)
        with pytest.raises(ValueError, match="must be a non-negative integer for type 'step'"):
            RangeQuery(type="step", end=invalid_val)

    def test_timestamp_bounds(self):
        rq = RangeQuery(type="timestamp", start=1700000000000, end=1700000001000)
        assert rq.start == 1700000000000.0

        with pytest.raises(ValueError, match="must be a non-negative integer for type 'timestamp'"):
            RangeQuery(type="timestamp", start=-100)
        with pytest.raises(ValueError, match="must be a non-negative integer for type 'timestamp'"):
            RangeQuery(type="timestamp", end=10.5)

    def test_custom_bounds_allow_negative_and_floats(self):
        rq = RangeQuery(type="custom", start=-10.5, end=1.25e-3)
        assert rq.start == -10.5
        assert rq.end == 1.25e-3

    def test_custom_bounds_start_greater_than_end_raises(self):
        with pytest.raises(ValueError, match="start must be <= end"):
            RangeQuery(type="custom", start=10.0, end=5.0)

    def test_head_and_tail_mutually_exclusive(self):
        with pytest.raises(ValueError, match="head and tail are mutually exclusive"):
            RangeQuery(type="step", head=10, tail=10)

    def test_last_and_start_end_mutually_exclusive(self):
        with pytest.raises(ValueError, match="last is mutually exclusive with start/end"):
            RangeQuery(type="step", last=1000, start=10)
        with pytest.raises(ValueError, match="last is mutually exclusive with start/end"):
            RangeQuery(type="step", last=1000, end=20)
