import enum
from datetime import datetime
from functools import wraps
from typing import Any, Callable, Dict, Optional, get_args

import click
import nanoid
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import (
    _VALID_PAGE_SIZES,
    ApiColumnClassLiteral,
    ApiColumnDataTypeLiteral,
    ApiMetricLogLevelLiteral,
    ApiResponseType,
    ApiVisibilityLiteral,
)


class _SaveFormatEnum(enum.Enum):
    JSON = "json"


PAGE_SIZE_TYPE = click.Choice([str(s) for s in _VALID_PAGE_SIZES])
COLUMN_CLASS_TYPE = click.Choice(list(get_args(ApiColumnClassLiteral)), case_sensitive=False)
COLUMN_DATA_TYPE = click.Choice(list(get_args(ApiColumnDataTypeLiteral)), case_sensitive=False)
VISIBILITY_TYPE = click.Choice(list(get_args(ApiVisibilityLiteral)), case_sensitive=False)
METRIC_LOG_LEVEL_TYPE = click.Choice(list(get_args(ApiMetricLogLevelLiteral)), case_sensitive=False)


def with_custom_host(func: Callable) -> Callable:
    """
    Add common SwanLab API host/auth options to a CLI command.

    The wrapped command receives an `api` keyword argument. When no option is
    provided, the default local login settings are used.
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
    @wraps(func)
    def wrapper(*args, host: Optional[str] = None, api_key: Optional[str] = None, **kwargs):
        if host is None and api_key is None:
            api = Api()
        else:
            api = Api(host=host, api_key=api_key)
        return func(*args, api=api, **kwargs)

    return wrapper


def format_output(resp: ApiResponseType, fmt: _SaveFormatEnum = _SaveFormatEnum.JSON) -> Dict[str, Any]:
    payload = resp.json()
    if fmt == _SaveFormatEnum.JSON:
        click.echo(orjson.dumps(payload, option=orjson.OPT_INDENT_2).decode())
    return payload


def save_output(content: bytes, name: Optional[str] = None, fmt: _SaveFormatEnum = _SaveFormatEnum.JSON) -> None:
    if name and name != ".":
        ext = name.rsplit(".", 1)[-1].lower() if "." in name else None
        if ext and ext not in {f.value for f in _SaveFormatEnum}:
            click.echo(f"Warning: unsupported file extension .{ext}, skipped saving.")
            return
        filename = name
    else:
        filename = f"swanlab-{datetime.now().strftime('%Y%m%d_%H%M%S')}-{nanoid.generate(size=4)}.{fmt.value}"
    with open(filename, "wb") as f:
        f.write(content)
    click.echo(f"Saved to {filename}")


def parse_keys(keys: str) -> list[str]:
    """Parse comma-separated keys string into a list, raising click.BadParameter on empty result."""
    key_list = [k.strip() for k in keys.split(",") if k.strip()]
    if not key_list:
        raise click.BadParameter("No valid keys provided. Expected comma-separated keys, e.g. 'loss,acc'.")
    return key_list
