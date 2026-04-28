import click
import orjson

from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import PAGE_SIZE_TYPE, format_output, save_output, with_custom_host


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
def list_experiments(page_num: int, page_size: str, project_path: str, fetch_all: bool, name, api):
    """List experiments under a project."""
    resp = ApiResponseType(ok=True, data=api.runs_get(path=project_path, page=page_num, size=int(page_size), all=fetch_all))
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)


@experiment_cli.command("columns")
@click.argument("path", required=True)
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
    "--class",
    "column_class",
    default="CUSTOM",
    type=str,
    help="Column class, such as CUSTOM or SYSTEM.",
)
@click.option(
    "--type",
    "column_type",
    default=None,
    type=str,
    help="Column type, such as IMAGE, FLOAT, or STRING.",
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
def list_experiment_columns(
    path: str,
    page_num: int,
    page_size: str,
    column_class: str,
    column_type: str,
    fetch_all: bool,
    name,
    api,
):
    """List columns under an experiment."""
    resp = ApiResponseType(
        ok=True,
        data=api.columns(
            path=path,
            page=page_num,
            size=int(page_size),
            column_class=column_class,
            column_type=column_type,
            all=fetch_all,
        ),
    )
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)
