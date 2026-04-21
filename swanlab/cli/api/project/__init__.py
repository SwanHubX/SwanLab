import click

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("project")
def project_cli():
    """Project management commands."""
    pass
