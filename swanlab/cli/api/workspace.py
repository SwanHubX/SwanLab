from typing import Optional

import click
import orjson

from swanlab.api.typings.common import ApiResponseType


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass


@workspace_cli.command("info")
@click.argument("username", required=True)
def get_info(username: Optional[str] = None):
    """Get workspace info."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.workspace(username)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())
