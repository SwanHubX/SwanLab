"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:09
@description: CLI Converter 模块：convert 命令，用于转换其他实验跟踪工具的日志
"""

import click

from swanlab.sdk.internal.pkg import safe


@click.command()
@click.option(
    "--type",
    "-t",
    "convert_type",
    default="tensorboard",
    type=click.Choice(["tensorboard", "wandb", "mlflow", "wandb-local"]),
    help="The type of the experiment tracking tool you want to convert from.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="SwanLab project name.",
)
@click.option(
    "--workspace",
    "-w",
    default=None,
    type=str,
    help="swanlab.init workspace parameter.",
)
@click.option(
    "--mode",
    default="online",
    type=click.Choice(["online", "local", "offline", "disabled"]),
    help="The mode of the swanlab run.",
)
@click.option(
    "-l",
    "--logdir",
    type=str,
    default=None,
    help="The directory where the swanlab log files are stored",
)
# Tensorboard options
@click.option("--tb-log-dir", type=str, default=None, help="The directory where the tensorboard log files are stored.")
@click.option("--tb-logdir", type=str, default=None, help="Deprecated: use --tb-log-dir instead.", hidden=True)
@click.option(
    "--tb-types",
    default="scalar",
    type=str,
    help="The types of the tensorboard log files to convert, support ['scalar', 'image', 'audio', 'text']. default 'scalar', split with ',' str.",
)
# wandb options
@click.option("--wb-project", type=str, default=None, help="The project name of the wandb runs.")
@click.option("--wb-entity", type=str, default=None, help="The entity name of the wandb runs.")
@click.option("--wb-runid", type=str, default=None, help="The run_id of the wandb run.")
# wandb-local options
@click.option("--wb-dir", type=str, default="./wandb", help="The directory where the wandb local log files are stored.")
@click.option("--wb-run-dir", type=str, default=None, help="The run directory of the wandb local log files.")
# mlflow options
@click.option(
    "--mlflow-url",
    type=str,
    default="http://127.0.0.1:5000",
    help="The tracking url of the mlflow server (default: http://127.0.0.1:5000).",
)
@click.option(
    "--mlflow-exp", type=str, default=None, help="The experiment 'name' or 'id' of the mlflow runs (required)."
)
@click.option("--mlflow-runid", type=str, default=None, help="The run id of a specific mlflow run to convert.")
# resume option
@click.option(
    "--resume",
    is_flag=True,
    default=False,
    help="Resume mode: sync source run_id to SwanLab (requires --wb-runid).",
)
def convert(
    convert_type: str,
    project: str,
    mode: str,
    workspace: str,
    logdir: str,
    tb_logdir: str,
    tb_log_dir: str,
    tb_types: str,
    wb_project: str,
    wb_entity: str,
    wb_runid: str,
    wb_dir: str,
    wb_run_dir: str,
    mlflow_url: str,
    mlflow_exp: str,
    mlflow_runid: str,
    resume: bool,
):
    """Convert the log files of other experiment tracking tools to SwanLab."""
    if tb_logdir is not None:
        click.echo("Warning: The option `--tb-logdir` is deprecated, use `--tb-log-dir` instead.")
        tb_log_dir = tb_logdir

    if resume and not wb_runid:
        raise click.UsageError("--resume requires --wb-runid to specify a single run to resume.")

    if convert_type == "tensorboard":
        from swanlab.converter import TFBConverter

        with safe.block(message="TensorBoard conversion failed"):
            TFBConverter(project=project, workspace=workspace, mode=mode, log_dir=logdir, types=tb_types).run(
                convert_dir=tb_log_dir or ".", depth=3
            )

    elif convert_type == "wandb":
        from swanlab.converter import WandbConverter

        if not wb_project:
            raise click.UsageError("--wb-project is required when using wandb online converter")

        with safe.block(message="W&B conversion failed"):
            WandbConverter(project=project, workspace=workspace, mode=mode, log_dir=logdir, resume=resume).run(
                wb_project=wb_project, wb_entity=wb_entity, wb_run_id=wb_runid
            )

    elif convert_type == "wandb-local":
        from swanlab.converter import WandbLocalConverter

        with safe.block(message="W&B local conversion failed"):
            WandbLocalConverter(project=project, workspace=workspace, mode=mode, log_dir=logdir, resume=resume).run(
                root_wandb_dir=wb_dir, wandb_run_dir=wb_run_dir, wb_run_id=wb_runid
            )

    elif convert_type == "mlflow":
        from swanlab.converter import MLFlowConverter

        if not mlflow_exp:
            raise click.UsageError("--mlflow-exp is required when using mlflow converter.")

        with safe.block(message="MLflow conversion failed"):
            MLFlowConverter(project=project, workspace=workspace, mode=mode, log_dir=logdir).run(
                tracking_uri=mlflow_url, experiment=mlflow_exp, run_id=mlflow_runid
            )
