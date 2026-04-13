"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:09
@description: CLI Converter 模块：convert 命令，用于转换其他实验跟踪工具的日志
"""

import click


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
    default="cloud",
    type=click.Choice(["cloud", "local", "offline", "disabled"]),
    help="The mode of the swanlab run.",
)
@click.option(
    "--logdir",
    "-l",
    type=str,
    default=None,
    help="The directory where the swanlab log files are stored.",
)
# tensorboard options
@click.option("--tb-logdir", type=str, default=None, help="The directory where the tensorboard log files are stored.")
@click.option(
    "--tb-types",
    default=None,
    type=str,
    help="The types of the tensorboard log files to convert, default is all types.",
)
# wandb options
@click.option("--wb-project", type=str, default=None, help="The project name of the wandb runs.")
@click.option("--wb-entity", type=str, default=None, help="The entity name of the wandb runs.")
@click.option("--wb-runid", type=str, default=None, help="The run_id of the wandb run.")
# mlflow options
@click.option("--mlflow-uri", type=str, default=None, help="The tracking uri of the mlflow runs.")
@click.option("--mlflow-exp", type=str, default=None, help="The experiment name or id of the mlflow runs.")
# wandb-local options
@click.option("--wb-dir", type=str, default="./wandb", help="The directory where the wandb local log files are stored.")
@click.option("--wb-run-dir", type=str, default=None, help="The run directory of the wandb local log files.")
def convert(
    convert_type: str,
    project: str,
    mode: str,
    workspace: str,
    logdir: str,
    tb_logdir: str,
    tb_types: str,
    wb_project: str,
    wb_entity: str,
    wb_runid: str,
    mlflow_uri: str,
    mlflow_exp: str,
    wb_dir: str,
    wb_run_dir: str,
    **kwargs,
):
    """Convert the log files of other experiment tracking tools to SwanLab."""
    # TODO: 接入重构后的 converter 逻辑
    pass
