#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/13 01:35
@File: watch.py
@IDE: pycharm
@Description:
    watch命令
"""
from swanlab.env import get_swanlog_dir, ROOT, PORT, HOST
from swanlab.log import swanlog
from swanboard import SwanBoardRun
import click
import os


@click.command()
@click.argument(
    "path",
    envvar=ROOT,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    nargs=1,
    required=False,
)
@click.option(
    "--host",
    "-h",
    default=lambda: os.environ.get(HOST, "127.0.0.1"),
    type=str,
    nargs=1,
    help="The host of swanlab web, default by 127.0.0.1",
)
@click.option(
    "--port",
    "-p",
    default=lambda: os.environ.get(PORT, 5092),
    nargs=1,
    type=click.IntRange(1, 65535),
    help="The port of swanlab web, default by 5092",
)
@click.option(
    "--logdir",
    "-l",
    default=lambda: os.environ.get(ROOT, None),
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    help="Specify the folder to store Swanlog, which is by default the folder where Swanlab Watch is run."
         "The option will be deprecated in the future, you can just use `swanlab watch <LOG PATH>` to specify the path."
)
@click.option(
    "--log-level",
    default="info",
    nargs=1,
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    help="The level of log, default by info; You can choose one of [debug, info, warning, error, critical]",
)
def watch(path: str, host: str, port: int, logdir: str, log_level: str):
    """
    Run this commands to turn on the swanlab service.
    """
    swanlog.set_level(log_level)
    # ----- 校验path，path如果被输入，已经由上层校验已存在，可读，是一个文件夹 -----
    if logdir is not None:
        swanlog.warning(
            "The option `--logdir` will be deprecated in the future, "
            f"you can just use `swanlab watch {logdir}` to specify the path."
        )
    # logdir 覆盖 path，接下来统一处理path而不再管logdir
    path = logdir if logdir is not None else path
    if path is not None:
        path = os.path.abspath(path)
        os.environ[ROOT] = path
    try:
        path = get_swanlog_dir()
    except ValueError as e:
        return click.BadParameter(str(e))
    except NotADirectoryError:
        return click.BadParameter("SWANLAB_LOG_DIR must be a directory")
    except FileNotFoundError:
        return click.BadParameter(f"The log folder `{path}` was not found")
    # ----- 校验host和port -----
    try:
        SwanBoardRun.is_valid_port(port)
        SwanBoardRun.is_valid_ip(host)
    except ValueError as e:
        return click.BadParameter(str(e))
    # ---- 启动服务 ----
    SwanBoardRun.run(
        path=path,
        host=host,
        port=port,
    )
