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
from swanlab.api.self_hosted import SelfHosted
from swanlab.api.series import Series
from swanlab.api.typings.common import PaginatedQuery, RangeQuery
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
            "create_settings",
            lambda: SimpleNamespace(api_key=None, api_host="https://api.swanlab.cn", web_host="https://swanlab.cn"),
        )

        with pytest.raises(AuthenticationError, match="No API key"):
            Api._resolve_credentials(None, None)

    def test_init_reuses_client_login_response(self, monkeypatch):
        """Client 初始化已完成鉴权，Api 不应再发一次登录请求。"""
        from swanlab.sdk.internal.pkg import scope as scope_mod

        api_module = importlib.import_module("swanlab.api")

        def fake_client(api_key, base_url):
            scope_mod.set_context("login_resp", {"userInfo": {"username": "alice", "name": "Alice"}})
            return MagicMock()

        login_calls = []
        monkeypatch.setattr(api_module, "Client", fake_client)
        monkeypatch.setattr(
            Api,
            "_resolve_credentials",
            staticmethod(lambda api_key, host: ("key", "https://api.test", "https://web.test")),
        )
        monkeypatch.setattr(
            "swanlab.sdk.internal.pkg.client.bootstrap.login_by_api_key",
            lambda **kwargs: login_calls.append(kwargs) or {},
        )

        instance = Api()

        assert instance.username == "alice"
        assert login_calls == []

    def test_init_prefers_current_client_login_response(self, monkeypatch):
        """外层 Scope 的旧登录响应不应覆盖本次 Client 的登录结果。"""
        from swanlab.sdk.internal.pkg import scope as scope_mod

        api_module = importlib.import_module("swanlab.api")

        def fake_client(api_key, base_url):
            scope_mod.set_context("login_resp", {"userInfo": {"username": "alice", "name": "Alice"}})
            return MagicMock()

        monkeypatch.setattr(api_module, "Client", fake_client)
        monkeypatch.setattr(
            Api,
            "_resolve_credentials",
            staticmethod(lambda api_key, host: ("key", "https://api.test", "https://web.test")),
        )

        with scope_mod.Scope() as outer:
            outer.set("login_resp", {"userInfo": {"username": "stale", "name": "Stale"}})
            instance = Api()

        assert instance.username == "alice"

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
    def test_columns_accept_documented_column_types(self, ctx, api, column_type):
        ctx.client.get.return_value = _api_response(
            {"cuid": "run-cuid", "slug": "run1", "name": "test-run", "createdAt": "2024-08-01T00:00:00Z"}
        )
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
            _api_response(
                {"cuid": "run-cuid", "slug": "run-slug", "name": "test-run", "createdAt": "2024-08-01T00:00:00Z"}
            ),
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")

        column = exp.column("loss")

        assert column.run_id == "run-cuid"

    def test_metrics_uses_run_cuid_for_metric_payload(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
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
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
        ]
        assert payload["projectId"] == "project-cuid"
        assert payload["columns"] == [{"experimentId": "run-cuid", "key": "loss", "createdAt": 1722470400}]

    def test_medias_uses_run_cuid_for_metric_payload(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
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
        assert payload["columns"] == [{"experimentId": "run-cuid", "key": "image", "createdAt": 1722470400}]

    def test_logs_uses_run_cuid_for_log_params(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
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
        assert kwargs["params"]["createdAt"] == 1722470400


# ---------------------------------------------------------------------------
# Column / Columns — 校验 + 4xx
# ---------------------------------------------------------------------------
class TestColumnValidation:
    def test_columns_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid column_type"):
            Columns(ctx, path="user/proj/run1", query=PaginatedQuery(), column_type="INVALID", exp_created_at=1)  # type: ignore

    def test_columns_invalid_class_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid column_class"):
            Columns(ctx, path="user/proj/run1", query=PaginatedQuery(), column_class="INVALID", exp_created_at=1)  # type: ignore

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

        columns = Columns(ctx, path="user/proj/run-slug", query=PaginatedQuery(), exp_created_at=1722470400)

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

        col = Column(ctx, path="user/proj/run-slug", key="loss", exp_created_at=1722470400)

        assert col.name == "loss"
        assert col.run_id == "run-cuid"
        assert [call.args[0] for call in ctx.client.get.call_args_list] == [
            "/project/user/proj/runs/run-slug",
            "/experiment/run-cuid/column",
        ]

    def test_column_project_id_fetches_project_lazily(self, ctx):
        item = {"key": "loss", "name": "loss", "type": "FLOAT", "class": "CUSTOM"}
        col = Column(ctx, path="user/proj/run1", key="loss", data=cast(Any, item), exp_created_at=1722470400)
        ctx.client.get.return_value = _api_response({"cuid": "project-cuid"})

        assert col.project_id == "project-cuid"
        assert [call.args[0] for call in ctx.client.get.call_args_list] == ["/project/user/proj"]


# ---------------------------------------------------------------------------
# Metric / Metrics — 校验
# ---------------------------------------------------------------------------
class TestMetricValidation:
    def test_metric_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric_type"):
            Metric(ctx, project_id="p1", run_id="r1", key="loss", metric_type="INVALID", created_at=1)

    def test_metric_invalid_log_level_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric log level"):
            Metric(ctx, project_id="p1", run_id="r1", key="LOG", metric_type="LOG", log_level="VERBOSE", created_at=1)  # type: ignore

    def test_metric_scalar_no_key_raises(self, ctx):
        with pytest.raises(ValueError, match="key is required"):
            Metric(ctx, project_id="p1", run_id="r1", key="", metric_type="SCALAR", created_at=1)

    def test_metrics_empty_keys_raises(self, ctx):
        with pytest.raises(ValueError, match="non-empty"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=[], metric_type="SCALAR", created_at=1)

    def test_metrics_invalid_type_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid metric_type"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=["loss"], metric_type=cast(Any, "INVALID"), created_at=1)

    @pytest.mark.parametrize("keys", ["loss", [""], None])
    def test_metrics_invalid_keys_raises(self, ctx, keys):
        with pytest.raises(ValueError, match="keys must be a non-empty list"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=cast(List[str], keys), metric_type="SCALAR", created_at=1)


# ---------------------------------------------------------------------------
# createdAt — 类型与精度兼容 + payload 注入
# ---------------------------------------------------------------------------
class TestCreatedAt:
    def test_metrics_payload_contains_created_at(self, ctx):
        ctx.client.get.side_effect = [
            _api_response(
                {
                    "cuid": "run-cuid",
                    "slug": "run-slug",
                    "name": "test-run",
                    "project_id": "project-cuid",
                    "createdAt": "2024-08-01T00:00:00.000Z",
                }
            )
        ]
        ctx.client.post.return_value = _api_response([{"key": "loss", "metrics": []}, {"key": "loss", "latest": {}}])
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.metrics(keys=["loss"])

        scalar_calls = [c for c in ctx.client.post.call_args_list if c.args[0] == "/house/metrics/scalar"]
        assert len(scalar_calls) == 1
        assert scalar_calls[0].kwargs["data"]["columns"] == [
            {"experimentId": "run-cuid", "key": "loss", "createdAt": 1722470400}
        ]

    def test_logs_params_contain_created_at(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
                )
            if path == "/house/metrics/log":
                return _api_response({"logs": [{"message": "ready"}], "count": 1})
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.logs()

        _, kwargs = ctx.client.get.call_args_list[-1]
        assert kwargs["params"]["createdAt"] == 1722470400

    def test_summary_payload_contains_created_at(self, ctx):
        ctx.client.get.side_effect = [
            _api_response(
                {
                    "cuid": "run-cuid",
                    "slug": "run-slug",
                    "name": "test-run",
                    "project_id": "project-cuid",
                    "createdAt": "2024-08-01T00:00:00Z",
                }
            )
        ]
        ctx.client.post.return_value = _api_response({"run-cuid": {"loss": {"latest": {"value": 1.0}}}})
        exp = Experiment(ctx, path="user/proj/run-slug")

        exp.summary()

        payload = ctx.client.post.call_args_list[0].kwargs["data"]
        assert payload["experiments"] == [
            {
                "projectId": "project-cuid",
                "experimentId": "run-cuid",
                "createdAt": 1722470400,
            }
        ]

    def test_payload_raises_when_created_at_missing(self, ctx):
        ctx.client.get.side_effect = [
            _api_response({"cuid": "run-cuid", "slug": "run-slug", "name": "test-run", "project_id": "project-cuid"})
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")

        # createdAt 缺失时必须强制报错，禁止静默降级为无下界查询（慢查询）
        with pytest.raises(ValueError, match="Failed to parse createdAt of experiment 'user/proj/run-slug'"):
            exp.summary()

    def test_created_at_ts_property(self, ctx):
        ctx.client.get.side_effect = [
            _api_response(
                {
                    "cuid": "run-cuid",
                    "slug": "run-slug",
                    "name": "test-run",
                    "project_id": "project-cuid",
                    "createdAt": "2024-08-01T00:00:00Z",
                }
            )
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")

        assert exp.created_at_ts == 1722470400

    def test_created_at_ts_property_raises_when_experiment_missing(self, ctx_404):
        exp = Experiment(ctx_404, path="user/proj/run123")

        with pytest.raises(ValueError, match="Failed to load experiment 'user/proj/run123'.*404"):
            _ = exp.created_at_ts

    def test_series_requests_contain_created_at(self, ctx):
        ctx.client.get.side_effect = [
            _api_response(
                {
                    "cuid": "run-cuid",
                    "slug": "run-slug",
                    "name": "test-run",
                    "project_id": "project-cuid",
                    "createdAt": "2024-08-01T00:00:00Z",
                }
            )
        ]
        ctx.client.post.return_value = _api_response({"keys": ["loss"], "hasMore": False, "nextCursor": ""})
        exp = Experiment(ctx, path="user/proj/run-slug")

        list(exp.series())

        keys_payload = ctx.client.post.call_args_list[0].kwargs["data"]
        assert keys_payload["experiments"] == [
            {"projectId": "project-cuid", "experimentId": "run-cuid", "createdAt": 1722470400}
        ]

    def test_series_raises_when_created_at_missing(self, ctx):
        ctx.client.post.return_value = _api_response({"keys": ["loss"], "hasMore": False, "nextCursor": ""})

        # experiment ref 缺失 createdAt 时必须强制报错，禁止静默降级为无下界查询（慢查询）
        with pytest.raises(ValueError, match="createdAt"):
            Series(
                ctx.client,
                experiments=[{"projectId": "project-cuid", "experimentId": "run-cuid"}],
            )

    def test_summary_raises_load_error_when_experiment_missing(self, ctx_404):
        exp = Experiment(ctx_404, path="user/proj/run123")

        # 实验不存在时必须抛出携带路径与原始服务端错误的异常，而非空字符串时间戳错误
        with pytest.raises(ValueError, match="Failed to load experiment 'user/proj/run123'.*404"):
            exp.summary()

    def test_metrics_raises_load_error_when_experiment_missing(self, ctx_404):
        exp = Experiment(ctx_404, path="user/proj/run123")

        with pytest.raises(ValueError, match="Failed to load experiment 'user/proj/run123'.*404"):
            exp.metrics(["loss"])


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


# ---------------------------------------------------------------------------
# Metrics 自定义 X 轴单测
# ---------------------------------------------------------------------------
class TestMetricsWithCustomXAxis:
    def test_metrics_sampled_with_custom_x_axis(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
                )
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        ctx.client.post.side_effect = [
            _api_response([{"metrics": [{"step": 1, "value": 0.5}], "key": "loss"}]),
            _api_response([{"min": {"value": 0.5}}]),
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")
        res = exp.metrics(["loss"], x_axis="epoch")

        payload = ctx.client.post.call_args_list[0].kwargs["data"]
        assert payload["xKey"] == "epoch"
        assert payload["xType"] == "step"

        loss_item = res["list"][0]
        assert loss_item["xAxis"] == {"type": "CUSTOM", "key": "epoch"}

    def test_metrics_with_x_axis_in_keys_identity(self, ctx):
        """当 x_axis='epoch' 且 'epoch' 在 keys 中，epoch 序列的 value == index。"""

        def get(path, **kwargs):
            return _api_response(
                {
                    "cuid": "run-cuid",
                    "slug": "run-slug",
                    "name": "test-run",
                    "project_id": "project-cuid",
                    "createdAt": "2024-08-01T00:00:00Z",
                }
            )

        ctx.client.get.side_effect = get
        ctx.client.post.side_effect = [
            _api_response(
                [
                    {"metrics": [{"index": 10.0, "data": 10.0}], "key": "epoch"},
                    {"metrics": [{"index": 10.0, "data": 0.5}], "key": "loss"},
                ]
            ),
            _api_response([{"min": {"value": 10.0}}, {"min": {"value": 0.5}}]),
        ]
        exp = Experiment(ctx, path="user/proj/run-slug")
        res = exp.metrics(["epoch", "loss"], x_axis="epoch")

        assert len(res["list"]) == 2
        epoch_series = res["list"][0]
        loss_series = res["list"][1]

        assert epoch_series["key"] == "epoch"
        for pt in epoch_series["metrics"]:
            assert pt["data"] == pt["index"]

        assert loss_series["key"] == "loss"
        assert loss_series["xAxis"] == {"type": "CUSTOM", "key": "epoch"}

    def test_metrics_rejects_time_axis_in_csv_mode(self, ctx):
        ctx.client.get.return_value = _api_response(
            {"cuid": "run-cuid", "project_id": "p-1", "createdAt": "2024-08-01T00:00:00Z"}
        )
        exp = Experiment(ctx, path="user/proj/run-slug")
        with pytest.raises(ValueError, match="is not supported in CSV mode"):
            exp.metrics(["loss"], all=True, x_axis="time")

    def test_metrics_rejects_custom_range_without_custom_x_axis(self, ctx):
        ctx.client.get.return_value = _api_response(
            {"cuid": "run-cuid", "project_id": "p-1", "createdAt": "2024-08-01T00:00:00Z"}
        )
        exp = Experiment(ctx, path="user/proj/run-slug")
        with pytest.raises(ValueError, match="requires a custom x axis"):
            exp.metrics(["loss"], range_query={"type": "custom", "start": 0}, x_axis="step")

    def test_metrics_retains_root_ids(self, ctx):
        ctx.client.post.side_effect = [
            _api_response([{"metrics": [{"step": 1, "value": 0.5}], "key": "loss"}]),
            _api_response([{"min": {"value": 0.5}}]),
        ]
        metrics = Metrics(
            ctx,
            project_id="p-1",
            run_id="r-1",
            keys=["loss"],
            metric_type="SCALAR",
            root_pro_id="root-p",
            root_exp_id="root-r",
            created_at=1700000000,
        )
        res = list(metrics)
        assert len(res) == 1
        assert res[0].json()["rootProId"] == "root-p"
        assert res[0].json()["rootExpId"] == "root-r"


# ---------------------------------------------------------------------------
# Histogram — OpenAPI 直方图读取
# ---------------------------------------------------------------------------
class TestHistogramMetrics:
    @staticmethod
    def _hist_entry(key: str, step: int = 1):
        return {
            "experimentId": "r-1",
            "key": key,
            "metrics": [{"step": step, "timestamp": 1700000000000, "values": [1, 2], "bins": [0.0, 1.0, 2.0]}],
        }

    def test_histogram_batch_aligns_keys_and_fills_missing(self, ctx):
        """House 按 columns 顺序返回、无数据列缺席：SDK 按 key 对齐并补空结果。"""
        # 请求 ["a", "b", "c"]，响应缺席 "b"
        ctx.client.post.return_value = _api_response([self._hist_entry("a"), self._hist_entry("c", step=9)])
        metrics = Metrics(
            ctx, project_id="p-1", run_id="r-1", keys=["a", "b", "c"], metric_type="HISTOGRAM", created_at=1700000000
        )

        res = [m.json() for m in metrics]

        assert [r["key"] for r in res] == ["a", "b", "c"]  # 顺序与请求一致
        assert res[0]["metrics"][0]["values"] == [1, 2]
        assert res[1]["metrics"] == []  # 缺席列补空
        assert res[2]["metrics"][0]["step"] == 9
        # json() 输出不含 scalar/media 专属字段
        for r in res:
            assert "logs" not in r and "count" not in r and "steps" not in r and "step" not in r

    def test_histogram_experiment_entrypoint_payload(self, ctx):
        def get(path, **kwargs):
            if path == "/project/user/proj/runs/run-slug":
                return _api_response(
                    {
                        "cuid": "run-cuid",
                        "slug": "run-slug",
                        "name": "test-run",
                        "project_id": "project-cuid",
                        "createdAt": "2024-08-01T00:00:00Z",
                    }
                )
            raise AssertionError(f"unexpected GET {path}")

        ctx.client.get.side_effect = get
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        exp = Experiment(ctx, path="user/proj/run-slug")

        # dict 入参走 RangeQuery 校验后同样映射到 range
        exp.histograms(["grad"], step_range={"start": 10, "end": 20})

        call = ctx.client.post.call_args_list[0]
        assert call.args[0] == "/house/metrics/histogram"
        payload = call.kwargs["data"]
        assert payload["projectId"] == "project-cuid"
        assert payload["columns"] == [{"experimentId": "run-cuid", "key": "grad", "createdAt": 1722470400}]
        assert payload["num"] == 1500  # sample 映射为服务端抽帧上限
        assert payload["range"] == [10.0, 20.0]

    def test_histogram_batches_keys_in_fours(self, ctx):
        """超过 4 个 key 自动分批：5 keys → 2 次请求（4 + 1），单请求体积被钳制。"""
        ctx.client.post.return_value = _api_response([self._hist_entry("k1")])
        metrics = Metrics(
            ctx,
            project_id="p-1",
            run_id="r-1",
            keys=["k1", "k2", "k3", "k4", "k5"],
            metric_type="HISTOGRAM",
            created_at=1700000000,
        )

        list(metrics)

        assert ctx.client.post.call_count == 2
        first, second = ctx.client.post.call_args_list
        assert len(first.kwargs["data"]["columns"]) == 4
        assert len(second.kwargs["data"]["columns"]) == 1

    def test_histogram_rejects_all_and_range_query(self, ctx):
        """无全量导出通道：all / range_query 必须显式报错，不允许静默降级。"""
        with pytest.raises(ValueError, match="all=True is not supported"):
            Metric(ctx, project_id="p1", run_id="r1", key="grad", metric_type="HISTOGRAM", created_at=1, all=True)

        with pytest.raises(ValueError, match="all=True is not supported"):
            Metrics(ctx, project_id="p1", run_id="r1", keys=["grad"], metric_type="HISTOGRAM", created_at=1, all=True)

        with pytest.raises(ValueError, match="range_query is only supported for SCALAR"):
            Metrics(
                ctx,
                project_id="p1",
                run_id="r1",
                keys=["grad"],
                metric_type="HISTOGRAM",
                created_at=1,
                range_query=RangeQuery(type="step"),
            )

    def test_histogram_single_metric_dispatch_and_json(self, ctx):
        """单 key Metric 经 dispatch 委托批量路径，json() 保留 metrics 且剥除无关字段。"""
        ctx.client.post.return_value = _api_response([self._hist_entry("grad", step=3)])
        metric = Metric(ctx, project_id="p-1", run_id="r-1", key="grad", metric_type="HISTOGRAM", created_at=1700000000)

        result = metric.json()

        assert ctx.client.post.call_args.args[0] == "/house/metrics/histogram"
        assert result["key"] == "grad"
        assert result["metric_type"] == "HISTOGRAM"
        assert result["metrics"][0]["step"] == 3
        assert result["metrics"][0]["bins"] == [0.0, 1.0, 2.0]

    def test_histogram_failed_request_degrades_to_empty(self, ctx):
        """旧版 House 无此端点（404 → data=None）时优雅降级为空结果。"""
        ctx.client.post.return_value = _api_response(None)
        metrics = Metrics(
            ctx, project_id="p-1", run_id="r-1", keys=["a", "b"], metric_type="HISTOGRAM", created_at=1700000000
        )

        res = [m.json() for m in metrics]

        assert len(res) == 2
        assert all(r["metrics"] == [] for r in res)

    def test_histogram_series_type_accepted(self, ctx):
        """Series 层接受 HISTOGRAM，keys/availability 路径自动指向 histogram 端点。"""
        ctx.client.post.return_value = _api_response({"keys": ["grad"], "hasMore": False, "nextCursor": ""})
        series = Series(
            ctx,
            experiments=[{"projectId": "p-1", "experimentId": "r-1", "createdAt": 1700000000}],
            metric_type="HISTOGRAM",
        )

        keys = list(series)

        assert ctx.client.post.call_args.args[0] == "/house/metrics/histogram/keys"
        assert [k.key for k in keys] == ["grad"]
        assert all(k.metric_type == "HISTOGRAM" for k in keys)
        assert series.json()["metricType"] == "HISTOGRAM"

    def test_histogram_sample_clamped_to_1500(self, ctx):
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        metrics = Metrics(
            ctx, project_id="p-1", run_id="r-1", keys=["grad"], metric_type="HISTOGRAM", sample=9999, created_at=1
        )

        list(metrics)

        assert ctx.client.post.call_args.kwargs["data"]["num"] == 1500

    def test_histogram_sample_floor_at_one(self, ctx):
        """sample<=0 会被 House 的 gte=1 拒绝，钳制为 1 避免请求失败后静默降级。"""
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        metrics = Metrics(
            ctx, project_id="p-1", run_id="r-1", keys=["grad"], metric_type="HISTOGRAM", sample=0, created_at=1
        )

        list(metrics)

        assert ctx.client.post.call_args.kwargs["data"]["num"] == 1

    def test_histogram_step_range_in_payload(self, ctx):
        """step_range 复用 RangeQuery，start/end 透传为 House 的 range=[min,max]。"""
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        metrics = Metrics(
            ctx,
            project_id="p-1",
            run_id="r-1",
            keys=["grad"],
            metric_type="HISTOGRAM",
            step_range=RangeQuery(start=100, end=200),
            created_at=1,
        )

        list(metrics)

        assert ctx.client.post.call_args.kwargs["data"]["range"] == [100.0, 200.0]

    def test_histogram_step_range_one_sided(self, ctx):
        """单侧缺省映射为 House 的 0（该侧不限制）。"""
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        metrics = Metrics(
            ctx,
            project_id="p-1",
            run_id="r-1",
            keys=["grad"],
            metric_type="HISTOGRAM",
            step_range=RangeQuery(end=200),
            created_at=1,
        )

        list(metrics)

        assert ctx.client.post.call_args.kwargs["data"]["range"] == [0.0, 200.0]

    def test_histogram_step_range_validation(self, ctx):
        """step_range 复用 RangeQuery 的数值校验，并拒绝直方图没有语义的字段。"""
        with pytest.raises(ValueError, match="only supports type='step'"):
            Metrics(
                ctx,
                project_id="p-1",
                run_id="r-1",
                keys=["grad"],
                metric_type="HISTOGRAM",
                step_range=RangeQuery(type="timestamp", start=1, end=2),
                created_at=1,
            )

        with pytest.raises(ValueError, match="only supports start/end"):
            Metrics(
                ctx,
                project_id="p-1",
                run_id="r-1",
                keys=["grad"],
                metric_type="HISTOGRAM",
                step_range=RangeQuery(start=1, end=2, head=10),
                created_at=1,
            )

        with pytest.raises(ValueError, match="step_range is only supported for HISTOGRAM"):
            Metrics(
                ctx,
                project_id="p-1",
                run_id="r-1",
                keys=["grad"],
                metric_type="SCALAR",
                step_range=RangeQuery(start=1, end=2),
                created_at=1,
            )

        # 数值规则由 RangeQuery 统一负责：负数、小数、逆序都被拒绝
        with pytest.raises(ValueError):
            RangeQuery(start=-1)
        with pytest.raises(ValueError):
            RangeQuery(start=1.5)
        with pytest.raises(ValueError):
            RangeQuery(start=10, end=5)

    def test_histogram_without_step_range_omits_range(self, ctx):
        ctx.client.post.return_value = _api_response([self._hist_entry("grad")])
        metrics = Metrics(ctx, project_id="p-1", run_id="r-1", keys=["grad"], metric_type="HISTOGRAM", created_at=1)

        list(metrics)

        assert "range" not in ctx.client.post.call_args.kwargs["data"]

    def test_histogram_concurrency_capped_at_two(self, ctx, monkeypatch):
        """9 个 key → 3 批，直方图并发单独钳制到 2（不跟随全局 MAX_CONCURRENT_COUNT）。"""
        from swanlab.api.metric import SafeThreadPoolExecutor as real_executor

        ctx.client.post.return_value = _api_response([])
        seen: List[int] = []

        def spy(*args, **kwargs):
            seen.append(kwargs["max_workers"])
            return real_executor(*args, **kwargs)

        monkeypatch.setattr("swanlab.api.metric.SafeThreadPoolExecutor", spy)
        metrics = Metrics(
            ctx,
            project_id="p-1",
            run_id="r-1",
            keys=[f"k{i}" for i in range(9)],
            metric_type="HISTOGRAM",
            created_at=1,
        )

        list(metrics)

        assert seen == [2]
