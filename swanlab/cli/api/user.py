from typing import Optional

import click
import orjson


@click.command("user")
@click.argument("username", required=False)
def get_user(username: Optional[str] = None):
    """Get user info."""
    from swanlab.api import Api

    api = Api()
    user = api.user(username)
    click.echo(orjson.dumps(user.to_dict(), option=orjson.OPT_INDENT_2).decode())
