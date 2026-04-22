import click

from swanlab.api import Api
from swanlab.cli.api.helper import with_save_option


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass


@workspace_cli.command("info")
@click.argument("username", required=True)
@with_save_option
def get_workspace(username: str):
    """Get Workspace info."""
    api = Api()
    return api.workspace(username).wrapper()
