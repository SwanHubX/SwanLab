import click

from swanlab.api import Api
from swanlab.cli.api.helper import with_save_option


@click.group("run")
def experiment_cli():
    """Experiment(Run) management commands."""
    pass


@experiment_cli.command("info")
@click.argument("path", required=True)
@with_save_option
def get_experiment(path: str):
    """Get Experiment(Run) info by path (username/project/run_id)."""
    api = Api()
    return api.run(path)
