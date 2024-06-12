#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 01:39:04
@File: swanlab\cli\main.py
@IDE: vscode
@Description:
    swanlab脚本命令的主入口
"""
from swanlab.package import get_package_version
import swanlab.cli.commands as C
import click


@click.group(invoke_without_command=True)
@click.version_option(get_package_version(), "--version", "-v", message="SwanLab %(version)s")
def cli():
    pass


# noinspection PyTypeChecker
cli.add_command(C.login)  # 登录
# noinspection PyTypeChecker
cli.add_command(C.logout)  # 登出

# # ---------------------------------- watch命令，开启本地后端服务 ----------------------------------
#
#
# @cli.command()
# # 控制服务发布的ip地址
# @click.option(
#     "--host",
#     "-h",
#     default=None,
#     type=str,
#     help="The host of swanlab web, default by 127.0.0.1",
#     callback=options_rule.is_valid_ip,
# )
# # 控制服务发布的端口，默认5092
# @click.option(
#     "--port",
#     "-p",
#     default=None,
#     type=int,
#     help="The port of swanlab web, default by 5092",
#     callback=options_rule.is_valid_port,
# )
# # 实验文件夹
# @click.option(
#     "--logdir",
#     "-l",
#     default=None,
#     type=str,
#     help="Specify the folder to store Swanlog, which is by default the folder where Swanlab Watch is run.",
#     callback=options_rule.is_valid_root_dir,
# )
# # 日志等级
# @click.option(
#     "--log-level",
#     default="info",
#     type=click.Choice(["debug", "info", "warning", "error", "critical"]),
#     help="The level of log, default by info; You can choose one of [debug, info, warning, error, critical]",
# )
# def watch(log_level: str, **kwargs):
#     """Run this commands to turn on the swanlab service."""
#     from swanboard import run
#
#     run(log_level)


# # ---------------------------------- 转换命令，用于转换其他实验跟踪工具 ----------------------------------
# @cli.command()
# @click.option(
#     "--type",
#     "-t",
#     default="tensorboard",
#     type=click.Choice(["tensorboard", "wandb"]),
#     help="The type of the experiment tracking tool you want to convert to.",
# )
# @click.option(
#     "--project",
#     "-p",
#     default=None,
#     type=str,
#     help="SwanLab project name.",
# )
# @click.option(
#     "--workspace",
#     "-w",
#     default=None,
#     type=str,
#     help="swanlab.init workspace parameter.",
# )
# @click.option(
#     "--cloud",
#     default=True,
#     type=bool,
#     help="swanlab.init cloud parameter.",
# )
# @click.option(
#     "--logdir",
#     "-l",
#     type=str,
#     help="The directory where the swanlab log files are stored.",
# )
# @click.option(
#     "--tb_logdir",
#     type=str,
#     help="The directory where the tensorboard log files are stored.",
# )
# @click.option(
#     "--wb-project",
#     type=str,
#     help="The project name of the wandb runs.",
# )
# @click.option(
#     "--wb-entity",
#     type=str,
#     help="The entity name of the wandb runs.",
# )
# @click.option(
#     "--wb-runid",
#     type=str,
#     help="The run_id of the wandb run.",
# )
# def convert(
#         type: str,
#         project: str,
#         cloud: bool,
#         workspace: str,
#         logdir: str,
#         tb_logdir: str,
#         wb_project: str,
#         wb_entity: str,
#         wb_runid: str,
#         **kwargs,
# ):
#     """Convert the log files of other experiment tracking tools to SwanLab."""
#     if type == "tensorboard":
#         from swanlab.converter import TFBConverter
#
#         tfb_converter = TFBConverter(
#             convert_dir=tb_logdir,
#             project=project,
#             workspace=workspace,
#             cloud=cloud,
#             logdir=logdir,
#         )
#         tfb_converter.run()
#
#     elif type == "wandb":
#         from swanlab.converter import WandbConverter
#
#         print(wb_project, wb_entity, wb_runid)
#
#         wb_converter = WandbConverter(
#             project=project,
#             workspace=workspace,
#             cloud=cloud,
#             logdir=logdir,
#         )
#         wb_converter.run(
#             wb_project=wb_project,
#             wb_entity=wb_entity,
#             wb_run_id=wb_runid,
#         )
#
#     else:
#         raise ValueError("The type of the experiment tracking tool you want to convert to is not supported.")


if __name__ == "__main__":
    cli()
