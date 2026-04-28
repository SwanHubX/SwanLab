import click
import orjson

from swanlab.cli.api.helper import format_output, save_output, with_custom_host


@click.group("project")
def project_cli():
    """Project management commands."""
    pass


@project_cli.command("info")
@click.argument("path", required=True)
@click.option(
    "--save",
    "-s",
    "name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_project(path: str, name, api):
    """Get project info by path (username/project)."""
    resp = api.project(path).wrapper()
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)
