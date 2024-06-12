#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 01:39:04
@File: swanlab\cli\main.py
@IDE: vscode
@Description:
    swanlab脚本命令的主入口
"""
import click
from swanlab.package import get_package_version, is_login
from swanlab.api.auth import terminal_login
from ..env import get_swanlab_folder
from ..utils import FONT
import shutil
from swanboard.run import options_rule


@click.group(invoke_without_command=True)
@click.version_option(get_package_version(), "--version", "-v", message="SwanLab %(version)s")
def cli():
    pass


# ---------------------------------- watch命令，开启本地后端服务 ----------------------------------


@cli.command()
# 控制服务发布的ip地址
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host of swanlab web, default by 127.0.0.1",
    callback=options_rule.is_valid_ip,
)
# 控制服务发布的端口，默认5092
@click.option(
    "--port",
    "-p",
    default=None,
    type=int,
    help="The port of swanlab web, default by 5092",
    callback=options_rule.is_valid_port,
)
# 实验文件夹
@click.option(
    "--logdir",
    "-l",
    default=None,
    type=str,
    help="Specify the folder to store Swanlog, which is by default the folder where Swanlab Watch is run.",
    callback=options_rule.is_valid_root_dir,
)
# 日志等级
@click.option(
    "--log-level",
    default="info",
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    help="The level of log, default by info; You can choose one of [debug, info, warning, error, critical]",
)
def watch(log_level: str, **kwargs):
    """Run this command to turn on the swanlab service."""
    from swanboard import run

    run(log_level)


# ---------------------------------- 登录命令，进行登录 ----------------------------------
@cli.command()
@click.option(
    "--relogin",
    "-r",
    is_flag=True,
    default=False,
    help="Relogin to the swanlab cloud, it will recover the token file.",
)
@click.option(
    "--api-key",
    "-k",
    default=None,
    type=str,
    help="If you prefer not to engage in command-line interaction to input the api key, "
         "this will allow automatic login.",
)
def login(api_key: str, relogin: bool, **kwargs):
    """Login to the swanlab cloud."""
    if not relogin and is_login():
        # 此时代表token已经获取，需要打印一条信息：已经登录
        command = FONT.bold("swanlab login --relogin")
        tip = FONT.swanlab("You are already logged in. Use `" + command + "` to force relogin.")
        return print(tip)
    # 进行登录，此时将直接覆盖本地token文件
    login_info = terminal_login(api_key)
    print(FONT.swanlab("Login successfully. Hi, " + FONT.bold(FONT.default(login_info.username))) + "!")


# ---------------------------------- 退出登录命令 ----------------------------------
@cli.command()
def logout(**kwargs):
    """Logout to the swanlab cloud."""
    command = FONT.bold("swanlab login")
    if is_login():
        # 如果已经是登录状态，那么则询问用户是否确认，如果确认则删除token文件夹
        confirm = input(FONT.swanlab("Are you sure you want to logout? (y/N): "))
        if confirm.lower() == "y":
            try:
                shutil.rmtree(get_swanlab_folder())
                return print(FONT.swanlab("Logout successfully. You can use `" + command + "` to login again."))
            except Exception as e:
                return print(FONT.swanlab("Logout failed. Please check if you have file operation permissions."))
        else:
            return print(FONT.swanlab("Logout canceled."))

    # 如果还未登录，则不做任何处理，并告知用户如何登录
    tip = FONT.swanlab("You are not logged in. If you want to login in, please use `" + command + "` to login.")
    return print(tip)


# ---------------------------------- 转换命令，用于转换其他实验跟踪工具 ----------------------------------
@cli.command()
@click.option(
    "--type",
    "-t",
    default="tensorboard",
    type=click.Choice(["tensorboard", "wandb"]),
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
        **kwargs,
):
    """Convert the log files of other experiment tracking tools to SwanLab."""
    if type == "tensorboard":
        from swanlab.converter import TFBConverter

        tfb_converter = TFBConverter(
            convert_dir=tb_logdir,
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        tfb_converter.run()

    elif type == "wandb":
        from swanlab.converter import WandbConverter

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

    else:
        raise ValueError("The type of the experiment tracking tool you want to convert to is not supported.")


if __name__ == "__main__":
    cli()
