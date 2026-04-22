import click

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("selfhosted")
def selfhosted_cli():
    """Self-hosted deployment management commands."""
    pass
