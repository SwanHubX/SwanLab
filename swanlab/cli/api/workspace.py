import click

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import api_command


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass


@workspace_cli.command("info")
@click.argument("username", required=True)
@api_command
def get_workspace(username: str, api: Api) -> ApiResponseType:
    """Get Workspace info."""
    return api.workspace(username).wrapper()
