"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:10
@description: CLI Sync 模块：同步本地数据到云端
"""

from pathlib import Path
from typing import Optional, Tuple

import click

from swanlab.sdk import Settings, pkg
from swanlab.sdk import sync as sync_cmd
from swanlab.sdk.cmd.sync import build_sync_range


@click.command()
@click.argument(
    "path",
    nargs=-1,
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    required=True,
)
@click.option(
    "--api-key",
    "-k",
    default=None,
    type=str,
    help="The API key to use for authentication. If not specified, uses the default API key.",
)
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host to sync the logs to. If not specified, uses the default host.",
)
@click.option(
    "--workspace",
    "-w",
    default=None,
    type=str,
    help="The workspace to sync the logs to. If not specified, uses the previous workspace.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="The project to sync the logs to. If not specified, uses the previous project.",
)
@click.option(
    "--id",
    "-i",
    default=None,
    type=str,
    help=(
        "The run id to sync the logs to. If it matches an existing experiment, the data is "
        "replayed into that experiment. Only valid with a single path."
    ),
)
@click.option(
    "--start-step",
    default=None,
    type=click.IntRange(min=0),
    help="Only sync scalar/media records whose step is greater than or equal to this value.",
)
@click.option(
    "--end-step",
    default=None,
    type=click.IntRange(min=0),
    help="Only sync scalar/media records whose step is less than or equal to this value.",
)
@click.option(
    "--start-time",
    default=None,
    type=str,
    help=(
        "Only sync records at or after this time. Accepts ISO 8601 (e.g. '2026-10-20 10:00:00') "
        "or a Unix timestamp in seconds/milliseconds. Times without offset use the system "
        "local timezone."
    ),
)
@click.option(
    "--end-time",
    default=None,
    type=str,
    help="Only sync records at or before this time. Same formats as --start-time.",
)
def sync(
    path: Tuple[str],
    api_key: Optional[str],
    workspace: Optional[str],
    project: Optional[str],
    host: Optional[str],
    id: Optional[str],
    start_step: Optional[int],
    end_step: Optional[int],
    start_time: Optional[str],
    end_time: Optional[str],
):
    """Synchronize local logs to the cloud.

    Passing any of --start-step/--start-time/--end-step/--end-time enables range
    sync, which only replays records inside the window. Sync always creates or
    resumes the target experiment; combine it with --id of an existing run to
    backfill that experiment. The local finish record decides the final state,
    and an existing experiment without one is left running. Times without an
    explicit offset use the system local timezone.
    """
    if id is not None and len(path) > 1:
        raise click.BadParameter("--id can only be used with a single path.", param_hint="--id")
    model = pkg.helper.strip_none(
        {
            "api_key": api_key,
            "api_host": host,
            "project": {"workspace": workspace, "name": project},
            "run": {"id": id},
        },
        strip_empty_str=True,
        strip_empty_dict=True,
    )
    settings = Settings.model_validate(model)
    try:
        build_sync_range(
            start_step=start_step,
            end_step=end_step,
            start_time=start_time,
            end_time=end_time,
        )
    except ValueError as e:
        raise click.BadParameter(str(e)) from e
    for p in path:
        sync_cmd(
            Path(p),
            settings=settings,
            start_step=start_step,
            end_step=end_step,
            start_time=start_time,
            end_time=end_time,
        )
