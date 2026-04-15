"""
@author: caddiesnew
@file: logout.py
@time: 2026/4/9 21:00
@description: swanlab.logout 方法，登出 SwanLab 平台
"""

from swanlab.sdk.cmd import utils
from swanlab.sdk.internal.pkg import console, nrc
from swanlab.sdk.internal.settings import settings
from swanlab.sdk.typings.cmd import LoginType

__all__ = ["logout_cli"]


def logout_cli(force: bool = False, save: LoginType = "root") -> bool:
    nrc_path = utils.get_nrc_path(save)
    # 1. 检查是否已登录
    if not nrc_path.exists():
        console.error("You are not logged in. Use `swanlab login` to login.")
        return False

    # 2. 交互式确认（非 force 模式下）
    if not force:
        if settings.interactive:
            console.info("Are you sure you want to logout? (y/N): ", end="")
            try:
                confirm = input()
            except (KeyboardInterrupt, EOFError):
                console.info("\nLogout canceled.")
                return False
            if confirm.lower() != "y":
                console.info("Logout canceled.")
                return False
        else:
            console.info("Use `swanlab.logout(force=True)` to skip confirmation in non-interactive mode.")
            return False
    nrc.remove(nrc_path)
    console.info("Logout successfully. You can use `swanlab login` to login again.")
    return True
