from swanboard import SwanBoardRun
import click
import os
from swanlab.env import ROOT


def click_error_wrapper(func):
    """装饰器：捕获 ValueError 并转为 click.BadParameter 异常"""

    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValueError as e:
            raise click.BadParameter(str(e))

    return wrapper


class BoardRunner(SwanBoardRun):
    """封装控制面板启动器"""

    @click_error_wrapper
    def is_valid_ip(self, ctx, param, ip: str) -> None:
        super().is_valid_ip(ip)

    @click_error_wrapper
    def is_valid_port(self, ctx, param, port: int) -> None:
        super().is_valid_port(port)

    def is_valid_root_dir(self, ctx, param, log_dir: str) -> str:
        """检测是否是合法的日志目录，保证其可读且存在

        Parameters
        ----------
        ctx : click.Context
            上下文
        param : click.Parameter
            参数
        log_dir : str
            带检测的日志目录
        """
        # 将日志目录注入环境变量，在这之前先转换为绝对路径

        if log_dir is None:
            return

        # 将传入的路径转换为绝对路径
        log_dir = os.path.abspath(log_dir)

        # 必须是一个绝对路径
        if not os.path.isabs(log_dir):
            raise ValueError("Log dir must be an absolute path: " + log_dir)
        # 路径必须存在
        if not os.path.isdir(log_dir):
            raise ValueError("Log dir is not a directory: " + log_dir)
        # 路径必须可读
        if not os.access(log_dir, os.R_OK):
            raise ValueError("Log dir is not readable: " + log_dir)

        os.environ[ROOT] = log_dir


board_runner = BoardRunner()
