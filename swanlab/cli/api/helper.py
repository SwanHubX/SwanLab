import json

import click

from swanlab.api.typings.common import ApiResponseType


def format_output(resp: ApiResponseType) -> None:
    """统一输出 ApiResponseType JSON。"""
    click.echo(json.dumps(resp.json()))
