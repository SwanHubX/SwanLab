import click

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output


@click.group("project")
def project_cli():
    """Project management commands."""
    pass


@project_cli.command("info")
@click.argument("path", required=True)
def get_project(path: str):
    """Get project info by path (username/project)."""
    api = Api()
    resp = api.project(path)
    format_output(resp)
