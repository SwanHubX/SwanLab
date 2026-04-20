from typing import Optional

import click
import orjson

from swanlab.sdk.internal.pkg import safe


@click.command("workspace")
@click.argument("username", required=False)
def get_workspace(username: Optional[str] = None):
    """Get workspace info."""
    from swanlab.api import Api

    api = Api()
    with safe.block(message="Failed to get worksapce info"):
        workspace = api.workspace(username)
        click.echo(orjson.dumps(workspace.to_dict(), option=orjson.OPT_INDENT_2).decode())
