from typing import Optional

import click
import orjson

from swanlab.api.typings.common import ApiResponseType


@click.command("user")
@click.argument("username", required=False)
def get_user(username: Optional[str] = None):
    """Get user info."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.user(username)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())
