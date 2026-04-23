import enum
from datetime import datetime
from typing import Optional

import click
import nanoid
import orjson

from swanlab.api.typings.common import ApiResponseType


class _SaveFormatEnum(enum.Enum):
    JSON = "json"


def format_output(resp: ApiResponseType, fmt: _SaveFormatEnum = _SaveFormatEnum.JSON) -> None:
    if fmt == _SaveFormatEnum.JSON:
        click.echo(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2).decode())


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
