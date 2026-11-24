"""
@author: cunyue
@file: greeting.py
@time: 2026/4/28 13:36
@description: 用于在终端打印run运行开始、结束的提示信息
"""

from typing import TYPE_CHECKING

from rich.text import Text

from swanlab.proto.swanlab.operation.v1.operation_pb2 import CLOUD_STATE_FAILED
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import console, fs
from swanlab.sdk.internal.pkg.helper.version import get_swanlab_version

if TYPE_CHECKING:
    from swanlab.sdk.internal.run import Run

__all__ = ["welcome", "goodbye", "print_cloud_link"]


def _print_version():
    version = get_swanlab_version()
    console.info("Tracking run with swanlab version", version)


def _print_save_dir(ctx: RunContext):
    # skip_store 仅 online 模式合法，此时本地不再创建 swanlog 日志目录，所有数据只传云端.
    if ctx.config.settings.core.skip_store:
        console.info("⚠️ ", Text("Local record store disabled; data will be uploaded to cloud only", "yellow"))
        return
    console.info("💾 Run data saved at", Text(str(ctx.config.run_dir), "magenta bold"))


def _print_online_tip(run: "Run"):
    url = run.url
    project_url = url.split("/runs/")[0]
    console.info("📁 View project at", Text(project_url, style=f"link {project_url} blue underline"))
    # 截断 run.name，超过 name_truncate_len 时，前后保留 20 个字符，中间显示 ...
    console.info(
        "🚀 View run",
        Text(fs.safe_truncate(run.name, max_length=45)[0], "yellow"),
        "at",
        Text(url, style=f"link {url} blue underline"),
    )


def print_cloud_link(run: "Run") -> None:
    """在云端路径延迟绑定成功后输出项目与实验链接。"""
    _print_online_tip(run)


def _print_recover_online_tip(ctx: RunContext) -> None:
    """在云端绑定尚未完成时提示后台恢复状态，避免输出未生效的链接。"""
    console.info(
        "⏳ SwanLab Cloud has not attached yet; recovering in background and saving data locally at",
        Text(str(ctx.config.run_dir), "magenta bold"),
    )


def _print_attach_failed_tip(ctx: RunContext, run: "Run", *, show_sync_hint: bool) -> None:
    """云端连接永久失败时说明原因；必要时提示如何上传本地记录。"""
    reason = getattr(run, "_cloud_error", None) or "cloud attach failed"
    console.warning("❌ SwanLab Cloud recovery failed:", Text(reason, "yellow"))
    if not show_sync_hint:
        return
    console.info(
        "💾 Data is saved at",
        Text(str(ctx.config.run_dir), "magenta bold"),
        "; after resolving the issue, upload it with",
        Text(f"`swanlab sync {ctx.config.run_dir}`", "bold"),
    )


def _print_local_tip(ctx: RunContext):
    console.info(
        "🌟 Run",
        Text(f"`swanlab watch {str(ctx.config.settings.log_dir)}`", "bold"),
        "to view SwanLab Experiment Dashboard",
    )


def _print_offline_tip(ctx: RunContext):
    console.info(
        "⚡️ Run",
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
        if run._path:
            _print_online_tip(run)
        elif run._cloud_state == CLOUD_STATE_FAILED:
            _print_attach_failed_tip(ctx, run, show_sync_hint=True)
        else:
            _print_recover_online_tip(ctx)
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
        if run._path:
            _print_online_tip(run)
        elif run._cloud_state == CLOUD_STATE_FAILED:
            _print_attach_failed_tip(ctx, run, show_sync_hint=False)
        else:
            _print_recover_online_tip(ctx)
    # 2. Local
    elif mode == "local":
        _print_local_tip(ctx)
    # 3. Offline
    elif mode == "offline":
        _print_offline_tip(ctx)
