#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/14 00:51
@File: __init__.py.py
@IDE: pycharm
@Description:
    转换命令，来自`swanlab.converter`，兼容其他可视化工具，转换为swanlab格式
    由于`swanlab.converter`可以被用户调用，所以转换代码作为一个模块封装，这里只是调用 `swanlab.converter` 的代码
"""

import click


@click.command()
@click.option(
    "--type",
    "-t",
    default="tensorboard",
    type=click.Choice(["tensorboard", "wandb", "mlflow"]),
    help="The type of the experiment tracking tool you want to convert to.",
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
    "--cloud",
    default=True,
    type=bool,
    help="swanlab.init cloud parameter.",
)
@click.option(
    "--logdir",
    "-l",
    type=str,
    help="The directory where the swanlab log files are stored.",
)
@click.option(
    "--tb_logdir",
    type=str,
    help="The directory where the tensorboard log files are stored.",
)
@click.option(
    "--wb-project",
    type=str,
    help="The project name of the wandb runs.",
)
@click.option(
    "--wb-entity",
    type=str,
    help="The entity name of the wandb runs.",
)
@click.option(
    "--wb-runid",
    type=str,
    help="The run_id of the wandb run.",
)
@click.option(
    "--mlflow-uri",
    type=str,
    help="The tracking uri of the mlflow runs.",
)
@click.option(
    "--mlflow-exp",
    type=str,
    help="The experiment name or id of the mlflow runs.",
)
    
def convert(
        type: str,
        project: str,
        cloud: bool,
        workspace: str,
        logdir: str,
        tb_logdir: str,
        wb_project: str,
        wb_entity: str,
        wb_runid: str,
        mlflow_uri: str,
        mlflow_exp: str,
        **kwargs,
):
    """Convert the log files of other experiment tracking tools to SwanLab."""
    if type == "tensorboard":
        from swanlab.converter.tfb import TFBConverter

        tfb_converter = TFBConverter(
            convert_dir=tb_logdir,
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        tfb_converter.run()

    elif type == "wandb":
        from swanlab.converter.wb import WandbConverter

        print(wb_project, wb_entity, wb_runid)

        wb_converter = WandbConverter(
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        wb_converter.run(
            wb_project=wb_project,
            wb_entity=wb_entity,
            wb_run_id=wb_runid,
        )

    elif type == "mlflow":
        from swanlab.converter.mlf import MLFLowConverter

        mlf_converter = MLFLowConverter(
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        
        mlf_converter.run(
            tracking_uri=mlflow_uri,
            experiment=mlflow_exp,
        )

    else:
        raise ValueError("The type of the experiment tracking tool you want to convert to is not supported.")
