import click

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("user")
def user_cli():
    """User management commands."""
    pass
