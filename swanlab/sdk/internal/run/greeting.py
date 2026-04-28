"""
@author: cunyue
@file: greeting.py
@time: 2026/4/28 13:36
@description: 用于在终端打印run运行开始、结束的提示信息
"""

from typing import TYPE_CHECKING

from rich.text import Text

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.helper.version import get_swanlab_version

if TYPE_CHECKING:
    from swanlab.sdk.internal.run import Run

__all__ = ["welcome", "goodbye"]


def _print_version():
    version = get_swanlab_version()
    console.info("Tracking run with swanlab version", version)


def _print_save_dir(ctx: RunContext):
    console.info("Run data will be saved locally in", Text(str(ctx.config.run_dir), "magenta bold"))


def _print_online_tip(run: "Run"):
    project_url = run.url.split("/runs/")[0]
    console.info(f"🏠 View project at {project_url}")
    console.info(f"🚀 View run at {run.url}")


def _print_local_tip(ctx: RunContext):
    console.info(
        "🌟 Run ",
        Text(f"`swanlab watch {str(ctx.config.settings.log_dir)}`", "bold"),
        "to view SwanLab Experiment Dashboard",
    )


def _print_offline_tip(ctx: RunContext):
    console.info(
        "☁️ Run ",
        Text(f"`swanlab sync {str(ctx.config.run_dir)}`", "bold"),
        "to sync this run",
    )


def welcome(ctx: RunContext, run: "Run"):
    """
    实验开始欢迎语，根据不同模式选择不同的欢迎语
    0. Disabled: 直接返回
    1. Online: 打印版本号、运行日志存储位置、实验名称、项目、实验的URL，并单独开启线程检查新版本
    2. Local: 打印版本号、运行日志存储位置、watch命令
    3. Offline: 答应版本号、运行日志存储位置、实验名称、sync命令
    """
    mode = ctx.config.settings.mode
    # 0. Disabled
    if mode == "disabled":
        return
    # 1. Online
    elif mode == "online":
        _print_version()
        _print_save_dir(ctx)
        if ctx.config.settings.experiment.name:
            console.info("Syncing run", Text(ctx.config.settings.experiment.name, "yellow"))
        _print_online_tip(run)
    # 2. Local
    elif mode == "local":
        _print_version()
        _print_save_dir(ctx)
        _print_local_tip(ctx)
    # 3. Offline
    elif mode == "offline":
        _print_version()
        _print_save_dir(ctx)
        _print_offline_tip(ctx)


def goodbye(ctx: RunContext, run: "Run"):
    mode = ctx.config.settings.mode
    # 0. Disabled
    if mode == "disabled":
        return
    # 1. Online
    elif mode == "online":
        _print_online_tip(run)
    # 2. Local
    elif mode == "local":
        _print_local_tip(ctx)
    # 3. Offline
    elif mode == "offline":
        _print_offline_tip(ctx)
