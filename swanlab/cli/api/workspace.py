import json
from typing import Optional

import click


@click.command("workspace")
@click.argument("username", required=False)
def get_workspace(username: Optional[str] = None):
    """Get workspace info."""
    from swanlab.api import Api

    api = Api()
    workspace = api.workspace(username)
    click.echo(json.dumps(workspace.to_dict(), indent=2, ensure_ascii=False))
