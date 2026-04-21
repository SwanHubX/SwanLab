import click

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("run")
def experiment_cli():
    """Experiment(Run) management commands."""
    pass
