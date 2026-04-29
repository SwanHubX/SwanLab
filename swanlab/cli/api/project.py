import click
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import PAGE_SIZE_TYPE, VISIBILITY_TYPE, format_output, save_output, with_custom_host


@click.group("project")
def project_cli():
    """Project management commands."""
    pass


@project_cli.command("info")
@click.argument("path", required=True)
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_project(path: str, save_name: str, api: Api):
    """Get project info by path (e.g. username/project_name)."""
    resp = api.project(path).wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@project_cli.command("list")
@click.option(
    "--page_num",
    "-n",
    default=1,
    type=click.IntRange(min=1),
    help="Page number.",
)
@click.option(
    "--page_size",
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
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_projects(page_num: int, page_size: str, workspace: str, fetch_all: bool, save_name: str, api: Api):
    """List projects under a workspace."""
    workspace = workspace or api.username
    resp = ApiResponseType(
        ok=True, data=api.projects(path=workspace, page=page_num, size=int(page_size), all=fetch_all)
    )
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@project_cli.command("create")
@click.option(
    "-n",
    "--name",
    required=True,
    type=str,
    help="Project name (1-100 chars, only 0-9a-zA-Z-_.+).",
)
@click.option(
    "-v",
    "--visibility",
    default="PRIVATE",
    type=VISIBILITY_TYPE,
    help="Project visibility, PUBLIC or PRIVATE. Default PRIVATE.",
)
@click.option(
    "-d",
    "--description",
    default=None,
    type=str,
    help="Project description.",
)
@click.option(
    "-w",
    "--workspace",
    default=None,
    type=str,
    help="Workspace username. Defaults to the current logged-in workspace.",
)
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    help="Save output as JSON to current directory.",
)
@with_custom_host
def create_project(name: str, visibility: str, description: str, workspace: str, save_name: str, api: Api):
    """Create a project in a workspace."""
    project = api.create_project(username=workspace, name=name, visibility=visibility.upper(), description=description)  # type: ignore
    if project is None:
        format_output(ApiResponseType(ok=False, errmsg="Failed to create project"))
        return
    resp = project.wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)
