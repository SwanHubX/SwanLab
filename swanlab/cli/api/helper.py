import json

import click

from swanlab.api.typings.common import ApiResponseType


def format_output(resp: ApiResponseType) -> None:
    """统一输出 ApiResponse JSON。"""
    click.echo(json.dumps(resp.to_json_dict()))
