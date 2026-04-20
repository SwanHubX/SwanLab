from typing import Optional

import click


@click.command("workspace")
@click.argument("username", required=False)
def get_workspace(username: Optional[str] = None):
    """Get workspace info."""
    raise NotImplementedError("cli api workspace")
