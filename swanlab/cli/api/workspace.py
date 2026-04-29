import click
import orjson

from swanlab.api import Api
from swanlab.cli.api.helper import format_output, save_output, with_custom_host


@click.group("workspace")
def workspace_cli():
    """Workspace management commands."""
    pass


@workspace_cli.command("info")
@click.argument("username", required=True)
@click.option(
    "--save",
    "-s",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_workspace(username: str, save_name: str, api: Api):
    """Get Workspace info."""
    resp = api.workspace(username).wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)
