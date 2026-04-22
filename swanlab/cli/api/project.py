import click

from swanlab.api import Api
from swanlab.cli.api.helper import with_save_option


@click.group("project")
def project_cli():
    """Project management commands."""
    pass


@project_cli.command("info")
@click.argument("path", required=True)
@with_save_option
def get_project(path: str):
    """Get project info by path (username/project)."""
    api = Api()
    return api.project(path)
