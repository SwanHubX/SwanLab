import click
import orjson

from swanlab.api.typings.common import ApiResponseType


def format_output(resp: ApiResponseType) -> None:
    """统一输出 ApiResponse JSON。"""
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())
