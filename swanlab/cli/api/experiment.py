import click
import orjson

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import format_output, save_output, with_custom_host


@click.group("experiment")
def experiment_cli():
    """Experiment(Run) management commands."""
    pass


@experiment_cli.command("info")
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
def get_experiment(path: str, name, api):
    """Get Experiment(Run) info by path (username/project/run_id)."""
    resp = api.run(path).wrapper()
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)


@experiment_cli.command("list")
@click.option(
    "--page_num",
    "--page-num",
    "-n",
    default=1,
    type=int,
    help="Page number.",
)
@click.option(
    "--page_size",
    "--page-size",
    "-s",
    default=20,
    type=int,
    help="Page size.",
)
@click.option(
    "--project_path",
    "--project-path",
    "-p",
    required=True,
    type=str,
    help="Project path, like username/project_name.",
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
def list_experiments(page_num: int, page_size: int, project_path: str, fetch_all: bool, name, api):
    """List experiments under a project."""
    resp = ApiResponseType(ok=True, data=api.runs_get(path=project_path, page=page_num, size=page_size, all=fetch_all))
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)
