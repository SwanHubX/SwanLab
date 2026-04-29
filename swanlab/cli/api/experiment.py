import click
import orjson

from swanlab.api import Api
from swanlab.api.typings.common import ApiResponseType
from swanlab.cli.api.helper import (
    COLUMN_CLASS_TYPE,
    COLUMN_DATA_TYPE,
    METRIC_LOG_LEVEL_TYPE,
    PAGE_SIZE_TYPE,
    format_output,
    parse_keys,
    save_output,
    with_custom_host,
)


@click.group("run")
def run_cli():
    """Experiment(Run) management commands."""
    pass


@run_cli.command("info")
@click.argument("path", required=True)
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_experiment(path: str, save_name: str, api: Api):
    """Get Experiment(Run) info by path (e.g. username/project_name/run_id)."""
    resp = api.run(path).wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("list")
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
    "--project_path",
    "-p",
    required=True,
    type=str,
    help="Project path (e.g. username/project_name).",
)
@click.option("--all", "fetch_all", is_flag=True, default=False, help="Fetch all pages.")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def list_experiments(page_num: int, page_size: str, project_path: str, fetch_all: bool, save_name: str, api: Api):
    """List experiments under a project by project path (e.g. username/project_name)."""
    resp = ApiResponseType(
        ok=True, data=api.runs_get(path=project_path, page=page_num, size=int(page_size), all=fetch_all)
    )
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("columns")
@click.argument("path", required=True)
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
    "--class",
    "column_class",
    default="CUSTOM",
    type=COLUMN_CLASS_TYPE,
    help="Column class, such as CUSTOM or SYSTEM.",
)
@click.option(
    "--type",
    "column_type",
    default=None,
    type=COLUMN_DATA_TYPE,
    help="Column type, such as IMAGE, FLOAT, or STRING.",
)
@click.option("--all", "fetch_all", is_flag=True, default=False, help="Fetch all pages.")
@click.option(
    "--save",
    "save_name",
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
    save_name: str,
    api: Api,
):
    """List columns under an experiment by path (e.g. username/project_name/run_id)."""
    resp = ApiResponseType(
        ok=True,
        data=api.columns(
            path=path,
            page=page_num,
            size=int(page_size),
            column_class=column_class.upper(),  # type: ignore
            column_type=column_type.upper() if column_type else None,  # type: ignore
            all=fetch_all,
        ),
    )
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("metrics")
@click.argument("path", required=True)
@click.option(
    "--keys",
    required=True,
    type=str,
    help="Comma-separated metric keys, e.g. 'loss,acc'.",
)
@click.option(
    "--sample",
    "-s",
    default=1500,
    type=click.IntRange(min=1),
    help="Sample size for scalar metrics. Default 1500.",
)
@click.option(
    "--ignore-timestamp",
    "ignore_timestamp",
    is_flag=True,
    default=False,
    help="Remove timestamp from metric data.",
)
@click.option("--all", "fetch_all", is_flag=True, default=False, help="Fetch all data (CSV export for scalars).")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_experiment_metrics(
    path: str,
    keys: str,
    sample: int,
    ignore_timestamp: bool,
    fetch_all: bool,
    save_name: str,
    api: Api,
):
    """Get scalar metrics of an experiment by path (e.g. username/project_name/run_id)."""
    key_list = parse_keys(keys)
    experiment = api.run(path)
    data = experiment.metrics(keys=key_list, sample=sample, ignore_timestamp=ignore_timestamp, all=fetch_all)
    payload = format_output(ApiResponseType(ok=True, data=data))
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("column")
@click.argument("path", required=True)
@click.option(
    "--key",
    required=True,
    type=str,
    help="Column key name.",
)
@click.option(
    "--class",
    "column_class",
    default="CUSTOM",
    type=COLUMN_CLASS_TYPE,
    help="Column class, such as CUSTOM or SYSTEM.",
)
@click.option(
    "--type",
    "column_type",
    default=None,
    type=COLUMN_DATA_TYPE,
    help="Column type, such as IMAGE, FLOAT, or STRING.",
)
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_experiment_column(
    path: str,
    key: str,
    column_class: str,
    column_type: str,
    save_name: str,
    api: Api,
):
    """Get column info of an experiment by path (e.g. username/project_name/run_id)."""
    col = api.column(
        path=path,
        key=key,
        column_class=column_class.upper(),  # type: ignore
        column_type=column_type.upper() if column_type else None,  # type: ignore
    )
    resp = col.wrapper()
    payload = format_output(resp)
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("medias")
@click.argument("path", required=True)
@click.option(
    "--keys",
    required=True,
    type=str,
    help="Comma-separated media keys, e.g. 'image,audio'.",
)
@click.option(
    "--step",
    "-s",
    default=0,
    type=int,
    help="Step number for media data. Default 0.",
)
@click.option("--all", "fetch_all", is_flag=True, default=False, help="Fetch all media steps.")
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_experiment_medias(
    path: str,
    keys: str,
    step: int,
    fetch_all: bool,
    save_name: str,
    api: Api,
):
    """Get media metrics of an experiment by path (e.g. username/project_name/run_id)."""
    key_list = parse_keys(keys)
    experiment = api.run(path)
    data = experiment.medias(keys=key_list, step=step, all=fetch_all)
    payload = format_output(ApiResponseType(ok=True, data=data))
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)


@run_cli.command("logs")
@click.argument("path", required=True)
@click.option(
    "--offset",
    "-o",
    default=0,
    type=int,
    help="Log offset (shard index). Default 0.",
)
@click.option(
    "--level",
    "-l",
    default="INFO",
    type=METRIC_LOG_LEVEL_TYPE,
    help="Log level: DEBUG, INFO, WARN, ERROR. Default INFO.",
)
@click.option(
    "--ignore-timestamp",
    "ignore_timestamp",
    is_flag=True,
    default=False,
    help="Remove timestamp from log data.",
)
@click.option(
    "--save",
    "save_name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
@with_custom_host
def get_experiment_logs(
    path: str,
    offset: int,
    level: str,
    ignore_timestamp: bool,
    save_name: str,
    api: Api,
):
    """Get console logs of an experiment by path (e.g. username/project_name/run_id)."""
    experiment = api.run(path)
    data = experiment.logs(offset=offset, level=level.upper(), ignore_timestamp=ignore_timestamp)  # type: ignore
    payload = format_output(ApiResponseType(ok=True, data=data))
    if payload["ok"] and save_name is not None:
        save_output(orjson.dumps(payload, option=orjson.OPT_INDENT_2), name=save_name)
