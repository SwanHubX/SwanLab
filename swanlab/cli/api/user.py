import click

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import api_command


@click.group("user")
def user_cli():
    """User management commands."""
    pass


@user_cli.command("info")
@api_command
def get_user(api: Api) -> ApiResponseType:
    """Get current User info"""
    return api.user().wrapper()
