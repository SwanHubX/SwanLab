from unittest.mock import MagicMock

import click
import orjson
import pytest
from click import BadParameter
from click.testing import CliRunner

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api import helper
from swanlab.cli.api.helper import validate_filter_query
from swanlab.exceptions import AuthenticationError


def test_inline_json():
    result = validate_filter_query('[{"key":"state","type":"STABLE","op":"EQ","value":["RUNNING"]}]')
    assert len(result) == 1
    assert result[0]["key"] == "state"


def test_from_file(tmp_path):
    p = tmp_path / "f.json"
    p.write_text('[{"key":"name","type":"STABLE","op":"CONTAIN","value":["test"]}]')
    result = validate_filter_query(str(p))
    assert len(result) == 1
    assert result[0]["op"] == "CONTAIN"


def test_empty_raises():
    with pytest.raises(BadParameter, match="must not be empty"):
        validate_filter_query("  ")


def test_invalid_json_raises():
    with pytest.raises(BadParameter, match="neither a valid file path nor valid JSON"):
        validate_filter_query("not json")


def test_non_array_raises():
    with pytest.raises(BadParameter, match="JSON array"):
        validate_filter_query('{"key": "state"}')


def test_api_command_success_passthrough(monkeypatch):
    """命令体正常返回时输出 ok=True 信封。"""
    monkeypatch.setattr(helper, "Api", MagicMock())

    @click.command()
    @click.argument("path")
    @helper.api_command
    def dummy(path, api):
        return ApiResponseType(ok=True, data={"path": path})

    result = CliRunner().invoke(dummy, ["x"])
    assert result.exit_code == 0
    assert '"ok": true' in result.output
    assert '"path": "x"' in result.output


def test_api_command_degrades_and_saves_on_error(monkeypatch, tmp_path):
    """命令体抛 ValueError（如实验不存在）应降级为 ok=False 信封，且 --save 仍保存完整 JSON。"""
    monkeypatch.setattr(helper, "Api", MagicMock())
    monkeypatch.chdir(tmp_path)

    @click.command()
    @click.argument("path")
    @helper.api_command
    def dummy(path, api):
        raise ValueError("Failed to load experiment 'user/proj/run': 404 Not Found")

    result = CliRunner().invoke(dummy, ["user/proj/run", "--save"])
    assert result.exit_code == 0  # 信封降级不算命令失败
    assert '"ok": false' in result.output
    assert "404 Not Found" in result.output
    saved = list(tmp_path.glob("swanlab-*.json"))
    assert len(saved) == 1
    assert orjson.loads(saved[0].read_bytes())["errmsg"].startswith("Failed to load experiment")


def test_api_command_converts_auth_error_to_click_exception(monkeypatch):
    """Api 构造阶段（未登录 / 认证失败）应转为单行 ClickException，而非 traceback。"""
    monkeypatch.setattr(helper, "Api", MagicMock(side_effect=AuthenticationError("No API key found")))

    @click.command()
    @helper.api_command
    def dummy():
        pass  # Api() 抛错时命令体不会执行

    result = CliRunner().invoke(dummy, [])
    assert result.exit_code == 1
    assert "Error: No API key found" in result.output
