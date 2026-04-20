import json
from typing import Optional

import click


@click.command("user")
@click.argument("username", required=False)
def get_user(username: Optional[str] = None):
    """Get user info."""
    from swanlab.api import Api

    api = Api()
    user = api.user(username)
    click.echo(json.dumps(user.to_dict(), indent=2, ensure_ascii=False))
