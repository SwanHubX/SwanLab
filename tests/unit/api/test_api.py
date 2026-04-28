"""
@author: caddiesnew
@time: 2026/4/27
@description: swanlab/api 实体类 4xx / 错误场景单测
"""

import importlib
from types import SimpleNamespace
from typing import Any, List, cast
from unittest.mock import MagicMock

import pytest
import requests

from swanlab.api import Api
from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.column import Column, Columns
from swanlab.api.experiment import Experiment, Experiments
from swanlab.api.metric import Metric, Metrics
from swanlab.api.project import Project
from swanlab.api.selfhosted import SelfHosted
from swanlab.api.typings.common import PaginatedQuery
from swanlab.api.typings.selfhosted import ApiSelfHostedInfoType
from swanlab.api.workspace import Workspace
from swanlab.exceptions import AuthenticationError

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MockResponse = MagicMock


def _api_response(data=None):
    """构造 Client.get/post 返回值。"""
    r = MagicMock()
    r.data = data
    r.raw = MagicMock()
    return r


@pytest.fixture
def ctx():
    client = MagicMock()
    return ApiClientContext(
        client=client,
        web_host="https://swanlab.cn",
        api_host="https://api.swanlab.cn",
        username="testuser",
        name="Test User",
    )


@pytest.fixture
def ctx_404(ctx):
    """Client 所有 HTTP 方法均抛出 HTTPError（模拟 4xx）。"""
    err = requests.exceptions.HTTPError("404 Not Found")
    ctx.client.get.side_effect = err
    ctx.client.post.side_effect = err
    ctx.client.put.side_effect = err
    ctx.client.delete.side_effect = err
    return ctx


@pytest.fixture
def api(ctx):
    instance = Api.__new__(Api)
    BaseEntity.__init__(instance, ctx)
    return instance


# ---------------------------------------------------------------------------
# Api 入口 — 参数校验
# ---------------------------------------------------------------------------
class TestApiEntryValidation:
    def test_missing_api_key_raises(self, monkeypatch):
        api_module = importlib.import_module("swanlab.api")
        monkeypatch.setattr(
            api_module,
            "global_settings",
            SimpleNamespace(api_key=None, api_host="https://api.swanlab.cn", web_host="https://swanlab.cn"),
        )

        with pytest.raises(AuthenticationError, match="No API key"):
            Api._resolve_credentials(None, None)

    @pytest.mark.parametrize("api_key", ["", "   "])
    def test_blank_api_key_raises(self, api_key):
        with pytest.raises(AuthenticationError, match="No API key"):
            Api._resolve_credentials(api_key, "https://api.swanlab.cn")

    def test_blank_host_raises(self):
        with pytest.raises(ValueError, match="Host cannot be empty"):
            Api._resolve_credentials("test-key", "   ")

    def test_projects_invalid_page_raises(self, api):
        with pytest.raises(ValueError, match="page must be >= 1"):
            api.projects("testuser", page=0)

    @pytest.mark.parametrize(
        ("method_name", "path"),
        [
            ("project", "testuser"),
            ("run", "testuser/project"),
            ("runs", "testuser/project/run1"),
            ("columns", "testuser/project"),
        ],
    )
    def test_factory_methods_reject_invalid_path_shapes(self, api, method_name, path):
        with pytest.raises(ValueError, match="path"):
            getattr(api, method_name)(path)

    def test_column_rejects_empty_key(self, api):
        with pytest.raises(ValueError, match="key"):
            api.column("testuser/project/run1", key="")

    @pytest.mark.parametrize("column_type", ["STRING", "IMAGE"])
    def test_columns_accept_documented_column_types(self, api, column_type):
        columns = api.columns("testuser/project/run1", column_type=column_type)
        assert isinstance(columns, Columns)


# ---------------------------------------------------------------------------
# SelfHosted — 权限拒绝
# ---------------------------------------------------------------------------
def _sh_info(**overrides) -> ApiSelfHostedInfoType:
    base = {"enabled": True, "expired": False, "root": True, "plan": "free", "seats": 10}
    base.update(overrides)
    return cast(ApiSelfHostedInfoType, base)


class TestSelfHostedPermission:
    def test_create_user_expired(self, ctx):
        sh = SelfHosted(ctx, data=_sh_info(expired=True))
        with pytest.raises(ValueError, match="expired"):
            sh.create_user("newuser", "pass123")

    def test_create_user_not_root(self, ctx):
        sh = SelfHosted(ctx, data=_sh_info(root=False))
        with pytest.raises(ValueError, match="root"):
            sh.create_user("newuser", "pass123")

    def test_get_users_not_root(self, ctx):
        sh = SelfHosted(ctx, data=_sh_info(root=False))
        with pytest.raises(ValueError, match="root"):
            list(sh.get_users())

    def test_create_user_4xx(self, ctx):
        sh = SelfHosted(ctx, data=_sh_info())
        ctx.client.post.side_effect = requests.exceptions.HTTPError("400 Bad Request")
        resp = sh.create_user("newuser", "pass123")
        assert not resp.ok

    @pytest.mark.parametrize(("username", "password"), [("", "pass123"), ("newuser", "")])
    def test_create_user_rejects_blank_credentials(self, ctx, username, password):
        sh = SelfHosted(ctx, data=_sh_info())
        with pytest.raises(ValueError, match="username|password"):
            sh.create_user(username, password)

    def test_get_users_4xx_yields_nothing(self, ctx):
        sh = SelfHosted(ctx, data=_sh_info())
        ctx.client.get.side_effect = requests.exceptions.HTTPError("500 Internal")
        result = list(sh.get_users())
        assert result == []


# ---------------------------------------------------------------------------
# Workspace — create_project 错误
# ---------------------------------------------------------------------------
class TestWorkspaceCreateProject:
    def test_invalid_name_raises_value_error(self, ctx):
        ws = Workspace(ctx, username="testuser")
        with pytest.raises(ValueError):
            ws.create_project("bad name!")

    def test_empty_name_raises_value_error(self, ctx):
        ws = Workspace(ctx, username="testuser")
        with pytest.raises(ValueError):
            ws.create_project("")

    def test_invalid_visibility_raises_value_error(self, ctx):
        ws = Workspace(ctx, username="testuser")
        with pytest.raises(ValueError):
            ws.create_project("valid-name", visibility=cast(Any, "SECRET"))
        ctx.client.post.assert_not_called()

    def test_api_error_returns_none(self, ctx):
        ws = Workspace(ctx, username="testuser")
        ctx.client.post.side_effect = requests.exceptions.HTTPError("500")
        result = ws.create_project("valid-name")
        assert result is None

    def test_success(self, ctx):
        ws = Workspace(ctx, username="testuser")
        ctx.client.post.return_value = _api_response(
            {"path": "testuser/valid-name", "name": "valid-name", "cuid": "cuid123"}
        )
        proj = ws.create_project("valid-name")
        assert proj is not None
        assert proj.name == "valid-name"


# ---------------------------------------------------------------------------
# Entity lazy-load 4xx — 确保不 crash，返回空默认值
# ---------------------------------------------------------------------------
class TestEntityLazyLoad4xx:
    def test_workspace_returns_empty(self, ctx_404):
        ws = Workspace(ctx_404, username="testuser")
        assert ws.name == ""
        assert ws.username == ""

    def test_project_returns_empty(self, ctx_404):
        proj = Project(ctx_404, path="user/proj")
        assert proj.name == ""
        assert proj.path == ""

    def test_experiment_returns_empty(self, ctx_404):
        exp = Experiment(ctx_404, path="user/proj/run123")
        assert exp.name == ""
        assert exp.state == ""

    def test_selfhosted_returns_defaults(self, ctx_404):
        sh = SelfHosted(ctx_404)
        assert sh.enabled is False
        assert sh.expired is False

    def test_project_delete_4xx_returns_false(self, ctx_404):
        proj = Project(ctx_404, path="user/proj")
        assert proj.delete() is False

    def test_experiment_delete_4xx_returns_false(self, ctx_404):
        exp = Experiment(ctx_404, path="user/proj/run123")
        assert exp.delete() is False


# ---------------------------------------------------------------------------
# Experiment — run slug / cuid 解析
# ---------------------------------------------------------------------------
class TestExperimentRunIdResolution:
    def test_run_id_resolves_slug_to_cuid(self, ctx):
        ctx.client.get.side_effect = [
            _api_response({"cuid": "run-cuid", "slug": "run-slug", "name": "test-run"}),
        ]

        exp = Experiment(ctx, path="user/proj/run-slug")

        assert exp.run_id == "run-cuid"
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
        ]

    def test_column_created_from_experiment_uses_run_cuid(self, ctx):
        ctx.client.get.side_effect = [
            _api_response({"cuid": "run-cuid", "slug": "run-slug", "name": "test-run"}),
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")

        column = exp.column("loss")

        assert column.run_id == "run-cuid"

    def test_metrics_uses_run_cuid_for_metric_payload(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {"cuid": "run-cuid", "slug": "run-slug", "name": "test-run", "project_id": "project-cuid"}
                )
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        ctx.client.post.side_effect = [
            _api_response([{"metrics": [{"step": 1, "value": 0.1}]}]),
            _api_response([{"min": {"value": 0.1}}]),
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.metrics(["loss"])

        payload = ctx.client.post.call_args_list[0].kwargs["data"]
        assert [call.args[0] for call in ctx.client.get.call_args_list] == ["/project/user/proj/runs/run-slug"]
        assert payload["projectId"] == "project-cuid"
        assert payload["columns"] == [{"experimentId": "run-cuid", "key": "loss"}]

    def test_medias_uses_run_cuid_for_metric_payload(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {"cuid": "run-cuid", "slug": "run-slug", "name": "test-run", "project_id": "project-cuid"}
                )
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        ctx.client.post.return_value = _api_response(
            {"steps": [1], "step": 1, "metrics": [{"key": "image", "data": [], "more": []}]}
        )
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.medias(["image"], step=1)

        payload = ctx.client.post.call_args_list[0].kwargs["data"]
        assert [call.args[0] for call in ctx.client.get.call_args_list] == ["/project/user/proj/runs/run-slug"]
        assert payload["projectId"] == "project-cuid"
        assert payload["columns"] == [{"experimentId": "run-cuid", "key": "image"}]

    def test_logs_uses_run_cuid_for_log_params(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {"cuid": "run-cuid", "slug": "run-slug", "name": "test-run", "project_id": "project-cuid"}
                )
            if path == "/house/metrics/log":
                return _api_response({"logs": [{"message": "ready"}], "count": 1})
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.logs()

        _, kwargs = ctx.client.get.call_args_list[-1]
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
            "/house/metrics/log",
        ]
        assert kwargs["params"]["projectId"] == "project-cuid"
        assert kwargs["params"]["experimentId"] == "run-cuid"


# ---------------------------------------------------------------------------
# Column / Columns — 校验 + 4xx
# ---------------------------------------------------------------------------
class TestColumnValidation:
    def test_columns_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid column_type"):
            Columns(ctx, path="user/proj/run1", query=PaginatedQuery(), column_type="INVALID")  # type: ignore

    def test_columns_invalid_class_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid column_class"):
            Columns(ctx, path="user/proj/run1", query=PaginatedQuery(), column_class="INVALID")  # type: ignore

    def test_iterated_columns_resolve_run_cuid_once_for_local_fields(self, ctx):
        ctx.client.get.side_effect = [
            _api_response({"cuid": "run-cuid", "slug": "run-slug"}),
            _api_response(
                {
                    "list": [
                        {"key": "loss", "name": "loss", "type": "FLOAT", "class": "CUSTOM"},
                        {"key": "acc", "name": "acc", "type": "FLOAT", "class": "CUSTOM"},
                    ],
                    "total": 2,
                    "pages": 1,
                }
            ),
        ]

        columns = Columns(ctx, path="user/proj/run-slug", query=PaginatedQuery())

        assert [column.name for column in columns] == ["loss", "acc"]
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
            "/experiment/run-cuid/column",
        ]

    def test_column_resolves_run_cuid_before_fetching_column_data(self, ctx):
        ctx.client.get.side_effect = [
            _api_response({"cuid": "run-cuid", "slug": "run-slug"}),
            _api_response(
                {
                    "list": [{"key": "loss", "name": "loss", "type": "FLOAT", "class": "CUSTOM"}],
                    "total": 1,
                    "pages": 1,
                }
            ),
        ]

        col = Column(ctx, path="user/proj/run-slug", key="loss")

        assert col.name == "loss"
        assert col.run_id == "run-cuid"
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
            "/experiment/run-cuid/column",
        ]

    def test_column_project_id_fetches_project_lazily(self, ctx):
        item = {"key": "loss", "name": "loss", "type": "FLOAT", "class": "CUSTOM"}
        col = Column(ctx, path="user/proj/run1", key="loss", data=cast(Any, item))
        ctx.client.get.return_value = _api_response({"cuid": "project-cuid"})

        assert col.project_id == "project-cuid"
        assert [call.args[0] for call in ctx.client.get.call_args_list] == ["/project/user/proj"]


# ---------------------------------------------------------------------------
# Metric / Metrics — 校验
# ---------------------------------------------------------------------------
class TestMetricValidation:
    def test_metric_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric_type"):
            Metric(ctx, project_id="p1", run_id="r1", key="loss", metric_type="INVALID")

    def test_metric_invalid_log_level_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric log level"):
            Metric(ctx, project_id="p1", run_id="r1", key="LOG", metric_type="LOG", log_level="VERBOSE")  # type: ignore

    def test_metric_scalar_no_key_raises(self, ctx):
        with pytest.raises(ValueError, match="key is required"):
            Metric(ctx, project_id="p1", run_id="r1", key="", metric_type="SCALAR")

    def test_metrics_empty_keys_raises(self, ctx):
        with pytest.raises(ValueError, match="non-empty"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=[], metric_type="SCALAR")

    def test_metrics_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric_type"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=["loss"], metric_type=cast(Any, "INVALID"))

    @pytest.mark.parametrize("keys", ["loss", [""], None])
    def test_metrics_invalid_keys_raises(self, ctx, keys):
        with pytest.raises(ValueError, match="keys must be a non-empty list"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=cast(List[str], keys), metric_type="SCALAR")


# ---------------------------------------------------------------------------
# Experiments POST 过滤 — 校验
# ---------------------------------------------------------------------------
class TestExperimentsFilterValidation:
    def test_invalid_filter_raises_on_iter(self, ctx):
        bad_filters = [{"key": "name"}]  # missing type, op, value
        exps = Experiments(ctx, path="user/proj", filters=bad_filters, mode="post")
        with pytest.raises(ValueError, match="Missing required"):
            list(exps)

    def test_non_list_filters_raise_on_iter(self, ctx):
        bad_filters = {"key": "name", "type": "STABLE", "op": "EQ", "value": ["test"]}
        exps = Experiments(ctx, path="user/proj", filters=cast(Any, bad_filters), mode="post")
        with pytest.raises(ValueError, match="filters must be a list"):
            list(exps)

    def test_invalid_group_raises_on_iter(self, ctx):
        bad_groups = [{"key": "cluster", "type": "INVALID"}]
        exps = Experiments(ctx, path="user/proj", groups=bad_groups, mode="post")
        with pytest.raises(ValueError, match="Invalid type"):
            list(exps)

    def test_invalid_sort_raises_on_iter(self, ctx):
        bad_sorts = [{"key": "name", "type": "STABLE", "order": "RANDOM"}]
        exps = Experiments(ctx, path="user/proj", sorts=bad_sorts, mode="post")
        with pytest.raises(ValueError, match="Invalid sort order"):
            list(exps)
