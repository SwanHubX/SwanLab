import functools
from datetime import datetime

import click
import nanoid
import orjson

from swanlab.api.typings.common import ApiResponseType


def _save_json(content: bytes) -> None:
    """将 JSON 内容保存到当前目录。"""
    filename = f"swanlab-{datetime.now().strftime('%Y%m%d_%H%M%S')}-{nanoid.generate(size=4)}.json"
    with open(filename, "wb") as f:
        f.write(content)
    click.echo(f"Saved to {filename}")


def format_output(resp: ApiResponseType, save: bool = False) -> None:
    """统一输出 ApiResponseType JSON，可选保存到文件。"""
    data = resp.json()
    click.echo(orjson.dumps(data, option=orjson.OPT_INDENT_2).decode())
    if save and resp.ok:
        _save_json(orjson.dumps(data, option=orjson.OPT_INDENT_2))


def with_save_option(f):
    """
    装饰器：为 CLI 命令添加 --save 选项并自动输出/保存响应。

    被装饰的函数应返回 ApiResponseType，装饰器负责 format_output 和可选的文件保存。
    """

    @click.option(
        "--save",
        "-s",
        is_flag=True,
        default=False,
        help="Save output as JSON to current directory.",
    )
    @functools.wraps(f)
    def wrapper(*args, save: bool, **kwargs):
        resp = f(*args, **kwargs)
        if resp is not None:
            format_output(resp, save=save)

    return wrapper
