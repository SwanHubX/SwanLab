import enum
from datetime import datetime
from functools import wraps
from typing import Any, Callable, Dict, Optional

import click
import nanoid
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType


class _SaveFormatEnum(enum.Enum):
    JSON = "json"


PAGE_SIZE_TYPE = click.Choice(["10", "12", "15", "20", "24", "27", "50", "100"])


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
    def wrapper(*args, host: Optional[str], api_key: Optional[str], **kwargs):
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
