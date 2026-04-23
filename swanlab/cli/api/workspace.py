import click
import orjson

from swanlab.api import Api
from swanlab.cli.api.helper import format_output, save_output


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass


@workspace_cli.command("info")
@click.argument("username", required=True)
@click.option(
    "--save",
    "-s",
    "name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
def get_workspace(username: str, name):
    """Get Workspace info."""
    api = Api()
    resp = api.workspace(username).wrapper()
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)
