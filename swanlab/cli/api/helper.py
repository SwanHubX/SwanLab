import enum
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, get_args

import click
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import (
    VALID_PAGE_SIZES,
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricKeyClassLiteral,
    ApiMetricKeyTypeLiteral,
    ApiMetricLogLevelLiteral,
    ApiResponseType,
    ApiVisibilityLiteral,
)
from swanlab.exceptions import ApiError, AuthenticationError
from swanlab.utils import generate_id


class _SaveFormatEnum(enum.Enum):
    JSON = "json"


PAGE_SIZE_TYPE = click.Choice([str(s) for s in VALID_PAGE_SIZES])
COLUMN_CLASS_TYPE = click.Choice(list(get_args(ApiColumnClassLiteral)), case_sensitive=False)
COLUMN_DATA_TYPE = click.Choice(list(get_args(ApiColumnDataTypeLiteral)), case_sensitive=False)
METRIC_KEY_TYPE = click.Choice(list(get_args(ApiMetricKeyTypeLiteral)), case_sensitive=False)
METRIC_KEY_CLASS = click.Choice(list(get_args(ApiMetricKeyClassLiteral)), case_sensitive=False)
VISIBILITY_TYPE = click.Choice(list(get_args(ApiVisibilityLiteral)), case_sensitive=False)
METRIC_LOG_LEVEL_TYPE = click.Choice(list(get_args(ApiMetricLogLevelLiteral)), case_sensitive=False)


def format_output(resp: ApiResponseType) -> bytes:
    content = orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2)
    click.echo(content)
    return content


def save_output(content: bytes, name: Optional[str] = None, fmt: _SaveFormatEnum = _SaveFormatEnum.JSON) -> None:
    if name and name != ".":
        ext = name.rsplit(".", 1)[-1].lower() if "." in name else None
        if ext and ext not in {f.value for f in _SaveFormatEnum}:
            click.echo(f"Warning: unsupported file extension .{ext}, skipped saving.")
            return
        filename = name
    else:
        filename = f"swanlab-{datetime.now().strftime('%Y%m%d_%H%M%S')}-{generate_id(length=4)}.{fmt.value}"
    with open(filename, "wb") as f:
        f.write(content)
    click.echo(f"Saved to {filename}")


def api_command(func: Callable) -> Callable:
    """
    SwanLab API CLI 命令统一装饰器。

    为CLI附加 ``--host`` / ``--api-key`` / ``--save`` 选项并注入已认证的 ``api`` 参数；

    - 正常返回： ``{"ok": true, "errmsg": "", "data": ...}``
    - 命令体抛出 ``ValueError``（实体不存在、createdAt 缺失等）：降级为 ``ok=False`` 响应并注入 ``errmsg``
    - 无论成功失败，``--save`` 均保存完整 JSON
    - Api 构造阶段（未登录 / 认证失败）不输出 JSON，通过单行 ClickException 提示认证失败
    """

    @click.option(
        "--host",
        "-h",
        default=None,
        type=str,
        help="The host of the SwanLab server.",
    )
    @click.option(
        "--api-key",
        "--api_key",
        "-k",
        "api_key",
        default=None,
        type=str,
        help="The API key to use for authentication.",
    )
    @click.option(
        "--save",
        "save_name",
        is_flag=False,
        flag_value=".",
        help="Save output as JSON to current directory.",
    )
    @wraps(func)
    def wrapper(
        *args, host: Optional[str] = None, api_key: Optional[str] = None, save_name: Optional[str] = None, **kwargs
    ):
        try:
            if host is None and api_key is None:
                api = Api()
            else:
                api = Api(host=host, api_key=api_key)
        except (ApiError, AuthenticationError) as e:
            raise click.ClickException(str(e)) from None
        try:
            resp = func(*args, api=api, **kwargs)
        except ValueError as e:
            resp = ApiResponseType(ok=False, errmsg=str(e))
        content = format_output(resp)
        if save_name is not None:
            save_output(content, name=save_name)
        return content

    return wrapper


def parse_keys(keys: str) -> list[str]:
    """Parse comma-separated keys string into a list, raising click.BadParameter on empty result."""
    key_list = [k.strip() for k in keys.split(",") if k.strip()]
    if not key_list:
        raise click.BadParameter("No valid keys provided. Expected comma-separated keys, e.g. 'loss,acc'.")
    return key_list


def validate_filter_query(query: str) -> List[Dict[str, Any]]:
    """
    Parse filter query from a file path or JSON string.

    If *query* points to an existing file, its contents are read and parsed as JSON.
    Otherwise it is treated as an inline JSON string.

    Returns a list of filter dicts (each must have key/type/op/value).
    """
    raw = query.strip()
    if not raw:
        raise click.BadParameter("filter_query must not be empty.")

    p = Path(raw)
    if p.is_file():
        try:
            data = orjson.loads(p.read_bytes())
        except (orjson.JSONDecodeError, OSError) as exc:
            raise click.BadParameter(f"Failed to read/parse filter file {raw!r}: {exc}")
    else:
        try:
            data = orjson.loads(raw)
        except orjson.JSONDecodeError as exc:
            raise click.BadParameter(f"filter_query is neither a valid file path nor valid JSON: {exc}")

    if not isinstance(data, list):
        raise click.BadParameter(f"filter_query must resolve to a JSON array, got {type(data).__name__}")

    return data
