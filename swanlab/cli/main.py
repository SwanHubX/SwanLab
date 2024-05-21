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
from .utils import is_valid_ip, is_valid_port, is_valid_root_dir, URL
from swanlab.package import version_limit, get_package_version, is_login
from swanlab.api.auth import terminal_login
from ..env import get_server_host, get_server_port, get_swanlog_dir
from ..db import connect
from ..utils import FONT
import time


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
    callback=is_valid_ip,
)
# 控制服务发布的端口，默认5092
@click.option(
    "--port",
    "-p",
    default=None,
    type=int,
    help="The port of swanlab web, default by 5092",
    callback=is_valid_port,
)
# 实验文件夹
@click.option(
    "--logdir",
    "-l",
    default=None,
    type=str,
    help="Specify the folder to store Swanlog, which is by default the folder where Swanlab Watch is run.",
    callback=is_valid_root_dir,
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
    # 本质上此模块用于注入环境变量并启动服务
    # 这里采用的是动态导入，因为路径的检查交由上面的is_valid_root_dir回调函数完成
    start = time.time()
    # 导入必要的模块
    from ..log import swanlog as swl
    from ..server import app
    import uvicorn

    log_dir = get_swanlog_dir()
    version_limit(log_dir, mode="watch")

    # debug一下当前日志文件夹的位置
    swl.debug("Try to explore the swanlab experiment logs in: " + FONT.bold(log_dir))
    try:
        connect()
    except FileNotFoundError:
        swl.error("Can not find the swanlab db in: " + FONT.bold(log_dir))
    # ---------------------------------- 日志等级处理 ----------------------------------
    swl.set_level(log_level)
    # ---------------------------------- 服务地址处理 ----------------------------------
    # 当前服务地址
    host = get_server_host()
    # 当前服务端口
    port = get_server_port()
    # 所有可用ip
    ipv4 = URL.get_all_ip()
    # ---------------------------------- 日志打印 ----------------------------------
    # 耗时
    take_time = int((time.time() - start) * 1000).__str__() + "ms\n\n"
    # 可用URL
    if URL.is_zero_ip(host):
        tip = "\n".join([URL(i, port).__str__() for i in ipv4])
    else:
        tip = URL(host, port).__str__()
    tip = tip + "\n" + URL.last_tip() + "\n"
    v = FONT.bold("v" + get_package_version())
    swl.info(f"SwanLab Experiment Dashboard " + v + " ready in " + FONT.bold(take_time) + tip)

    # ---------------------------------- 启动服务 ----------------------------------
    # 使用 uvicorn 启动 FastAPI 应用，关闭原生日志
    # 使用try except 捕获退出，涉及端口占用等
    try:
        uvicorn.run(app, host=host, port=port, log_level="critical")
    except SystemExit as e:
        code = e.code
        if code == 1:
            critical = "Error while attempting to bind on address ({}, {}): address already in use".format(host, port)
            swl.critical(critical)
        else:
            swl.critical("Unhandled Exit Code: {}".format(code))


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
