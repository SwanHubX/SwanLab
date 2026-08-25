from types import SimpleNamespace
from unittest.mock import MagicMock, call

import pytest

from swanlab.deprecated import project as deprecated_project_api
from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python.api import project as project_api


class _FakeApiErrorResponse:
    def __init__(self, status_code: int) -> None:
        self.status_code = status_code
        self.request = MagicMock()
        self.url = "https://api.test/project"


def _api_error(status_code: int, method: str = "POST") -> ApiError:
    return ApiError(
        _FakeApiErrorResponse(status_code),
        method=method,
        trace_id="trace-id",
        code="api-error",
        message="api error",
    )


def _project_data(username: str = "alice", name: str = "demo") -> dict:
    return {
        "cuid": "project-id",
        "name": name,
        "username": username,
        "path": f"/{username}/{name}",
        "visibility": "PRIVATE",
        "_count": {"experiments": 3, "contributors": 0, "collaborators": 0, "clones": 0},
    }


def _not_found_then_success(data: dict) -> MagicMock:
    """模拟项目详情接口：第一次 GET 404（项目不存在），第二次 GET 返回项目数据"""
    return MagicMock(side_effect=[_api_error(404, method="GET"), SimpleNamespace(data=data)])


def test_get_or_create_project_defaults_to_current_username(monkeypatch):
    post = MagicMock()
    get = _not_found_then_success(_project_data())
    username = MagicMock(return_value="alice")
    monkeypatch.setattr(project_api.client, "post", post)
    monkeypatch.setattr(project_api.client, "get", get)
    monkeypatch.setattr(project_api.client, "username", username)

    result = project_api.get_or_create_project(username=None, name="demo", public=False)

    assert result == _project_data()
    username.assert_called_once_with()
    post.assert_called_once_with(
        "/projects/alice",
        data={"name": "demo", "visibility": "PRIVATE"},
    )
    assert get.call_args_list == [call("/project/alice/demo", log_error=False), call("/project/alice/demo")]


def test_get_or_create_project_uses_explicit_username_without_current_username(monkeypatch):
    post = MagicMock()
    get = _not_found_then_success(_project_data(username="team", name="demo"))
    username = MagicMock(return_value="alice")
    monkeypatch.setattr(project_api.client, "post", post)
    monkeypatch.setattr(project_api.client, "get", get)
    monkeypatch.setattr(project_api.client, "username", username)

    result = project_api.get_or_create_project(username="team", name="demo", public=True)

    assert result == _project_data(username="team", name="demo")
    username.assert_not_called()
    post.assert_called_once_with(
        "/projects/team",
        data={"name": "demo", "visibility": "PUBLIC"},
    )
    assert get.call_args_list == [call("/project/team/demo", log_error=False), call("/project/team/demo")]


def test_get_or_create_project_existing_project_skips_create(monkeypatch):
    """项目已存在（GET 直接命中）时不应调用创建接口，避免只写权限用户触发 403"""
    post = MagicMock()
    get = MagicMock(return_value=SimpleNamespace(data=_project_data(username="team", name="demo")))
    monkeypatch.setattr(project_api.client, "post", post)
    monkeypatch.setattr(project_api.client, "get", get)

    result = project_api.get_or_create_project(username="team", name="demo", public=False)

    assert result == _project_data(username="team", name="demo")
    post.assert_not_called()
    get.assert_called_once_with("/project/team/demo", log_error=False)


def test_get_or_create_project_falls_back_to_old_helper_when_new_endpoint_fails(monkeypatch):
    post = MagicMock(side_effect=_api_error(404))
    get_or_create_old_project = MagicMock()
    get = _not_found_then_success(_project_data(username="team", name="demo"))
    monkeypatch.setattr(project_api.client, "post", post)
    monkeypatch.setattr(project_api, "get_or_create_old_project", get_or_create_old_project)
    monkeypatch.setattr(project_api.client, "get", get)

    result = project_api.get_or_create_project(username="team", name="demo", public=False)

    assert result == _project_data(username="team", name="demo")
    post.assert_called_once_with(
        "/projects/team",
        data={"name": "demo", "visibility": "PRIVATE"},
    )
    get_or_create_old_project.assert_called_once_with(
        data={"name": "demo", "visibility": "PRIVATE", "username": "team"}
    )
    assert get.call_args_list == [call("/project/team/demo", log_error=False), call("/project/team/demo")]


def test_get_or_create_project_raises_non_not_found_new_endpoint_api_error(monkeypatch):
    post = MagicMock(side_effect=_api_error(500))
    get_or_create_old_project = MagicMock()
    get = MagicMock(side_effect=_api_error(404, method="GET"))
    monkeypatch.setattr(project_api.client, "post", post)
    monkeypatch.setattr(project_api, "get_or_create_old_project", get_or_create_old_project)
    monkeypatch.setattr(project_api.client, "get", get)

    with pytest.raises(ApiError):
        project_api.get_or_create_project(username="team", name="demo", public=False)

    get_or_create_old_project.assert_not_called()
    get.assert_called_once_with("/project/team/demo", log_error=False)


def test_get_or_create_old_project_ignores_conflict(monkeypatch):
    post = MagicMock(side_effect=_api_error(409))
    monkeypatch.setattr(deprecated_project_api.client, "post", post)

    with pytest.warns(DeprecationWarning, match="legacy project"):
        deprecated_project_api.get_or_create_old_project(
            data={"name": "demo", "visibility": "PRIVATE", "username": "team"}
        )

    post.assert_called_once_with(
        "/project",
        data={"name": "demo", "visibility": "PRIVATE", "username": "team"},
        log_error=False,
    )
