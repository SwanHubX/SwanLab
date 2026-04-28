import click
import orjson

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import PAGE_SIZE_TYPE, format_output, save_output, with_custom_host


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
    payload = format_output(resp)
    if payload["ok"] and name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=name)


@project_cli.command("list")
@click.option(
    "--page_num",
    "--page-num",
    "-n",
    default=1,
    type=click.IntRange(min=1),
    help="Page number.",
)
@click.option(
    "--page_size",
    "--page-size",
    "-s",
    default="20",
    type=PAGE_SIZE_TYPE,
    help="Page size.",
)
@click.option(
    "--workspace",
    default=None,
    type=str,
    help="Workspace username. Defaults to the current logged-in workspace.",
)
@click.option("--all", "fetch_all", is_flag=True, default=False, help="Fetch all pages.")
@click.option(
    "--save",
    "name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_projects(page_num: int, page_size: str, workspace: str, fetch_all: bool, name, api):
    """List projects under a workspace."""
    workspace = workspace or api.username
    resp = ApiResponseType(ok=True, data=api.projects(path=workspace, page=page_num, size=int(page_size), all=fetch_all))
    payload = format_output(resp)
    if payload["ok"] and name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=name)
