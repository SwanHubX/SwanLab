"""
@author: caddiesnew
@time: 2026/4/27
@description: swanlab/api 校验函数单测
"""

from typing import cast

import pytest

from swanlab.api.selfhosted import SelfHosted
from swanlab.api.typings.common import PaginatedQuery
from swanlab.api.typings.selfhosted import ApiSelfHostedInfoType
from swanlab.api.utils import (
    validate_column_params,
    validate_filter,
    validate_group,
    validate_metric_log_level,
    validate_metric_type,
    validate_project_name,
    validate_sort,
)


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
        assert q.page == 1 and q.size == 20

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
