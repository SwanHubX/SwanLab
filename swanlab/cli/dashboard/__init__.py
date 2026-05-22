"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:09
@description: CLI Dashboard 模块：watch 命令，依赖于 swanboard
"""

import os
import socket
import sys

import click

from swanlab.sdk.internal import pkg


def _get_free_port(address: str = "0.0.0.0", default_port: int = 5092) -> int:
    """获取一个可用端口。默认返回 5092，如果被占用，返回一个随机可用端口。"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind((address, default_port))
        except OSError:
            pass
        else:
            return default_port
    sock = socket.socket()
    sock.bind((address, 0))
    _, port = sock.getsockname()
    sock.close()
    return port


@click.command()
@click.argument(
    "path",
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    nargs=1,
    required=False,
    default=None,
)
@click.option(
    "--host",
    "-h",
    default="127.0.0.1",
    type=str,
    nargs=1,
    help="The host of swanlab web, default by 127.0.0.1",
)
@click.option(
    "--port",
    "-p",
    default=None,
    nargs=1,
    type=click.IntRange(1, 65535),
    help="The port of swanlab web, default by 5092",
)
@click.option(
    "--log-dir",
    "-l",
    default=None,
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    help="Specify the folder to store Swanlog. Deprecated: use `swanlab watch <LOG PATH>` instead.",
)
@click.option(
    "--logdir",
    default=None,
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    help="Deprecated: use --log-dir instead.",
    hidden=True,
)
def watch(path: str, host: str, port: int, log_dir: str, logdir: str):
    """Run this command to turn on the SwanLab dashboard service."""
    # logdir 兼容旧参数（向后兼容）
    if logdir is not None:
        click.echo("Warning: The option `--logdir` is deprecated, use `--log-dir` instead.")
        log_dir = logdir

    # log-dir 覆盖 path（向后兼容）
    if log_dir is not None:
        click.echo(
            "Warning: The option `--log-dir` is deprecated, use `swanlab watch [PATH]` to specify the path instead."
        )
        path = log_dir

    if path is not None:
        path = os.path.abspath(path)

    if port is None:
        port = _get_free_port()

    try:
        # noinspection PyPackageRequirements
        from swanboard import SwanBoardRun
        from swanboard.utils import get_swanlog_dir
    except ModuleNotFoundError:
        raise
        click.echo("Please install the swanboard package: `pip install swanlab[dashboard]`")
        return sys.exit(1)
    # ----- 校验path，path如果被输入，已经由上层校验已存在，可读，是一个文件夹 -----
    if logdir is not None:
        pkg.console.warning(
            "The option `--logdir` will be deprecated in the future, "
            "you can just use `swanlab watch [PATH]` to specify the path."
        )
    # logdir 覆盖 path，接下来统一处理path而不再管logdir
    path = logdir if logdir is not None else path
    if path is not None:
        path = os.path.abspath(path)
        os.environ["SWANLAB_LOG_DIR"] = path
    # 为None时从环境变量中获取
    try:
        path = get_swanlog_dir()
        # 产品经理要求无论日志文件夹是否存在都不报错 🤡 ，给用户以开启web服务的“爽感”
        # if not os.path.exists(path):
        #     raise FileNotFoundError
    except ValueError as e:
        click.BadParameter(str(e))
        return sys.exit(3)
    except NotADirectoryError:
        click.BadParameter("SWANLAB_LOG_DIR must be a directory")
        return sys.exit(4)
    except FileNotFoundError:
        click.BadParameter(f"The log folder `{path}` was not found")
        return sys.exit(5)
    # ----- 校验host和port -----
    try:
        SwanBoardRun.is_valid_port(port)
        SwanBoardRun.is_valid_ip(host)
    except ValueError as e:
        click.BadParameter(str(e))
        return sys.exit(6)
    # ---- 启动服务 ----
    SwanBoardRun.run(
        path=path,
        host=host,
        port=port,
    )
