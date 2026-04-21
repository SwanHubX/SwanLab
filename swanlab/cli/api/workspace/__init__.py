import click

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass
