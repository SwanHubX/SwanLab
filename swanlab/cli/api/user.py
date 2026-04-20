from typing import Optional

import click


@click.command("user")
@click.argument("username", required=False)
def get_user(username: Optional[str] = None):
    """Get user info."""
    raise NotImplementedError("cli api user")
