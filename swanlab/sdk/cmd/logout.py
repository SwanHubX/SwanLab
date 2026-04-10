"""
@author: caddiesnew
@file: logout.py
@time: 2026/4/9 21:00
@description: swanlab.logout 方法，登出 SwanLab 平台
"""

from swanlab.sdk.cmd.helper import with_cmd_lock, without_run
from swanlab.sdk.internal import apikey
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.settings import settings


@with_cmd_lock
@without_run("logout")
def logout(force: bool = False) -> bool:
    """Logout from SwanLab Cloud.

    This function removes locally stored credentials and resets the runtime client.
    If no active login is found, the function will exit with an error.

    :param force: If True, skip the confirmation prompt and logout directly.
        Defaults to False.

    :return: True if logout was successful, False otherwise.

    :raises RuntimeError: If called while a run is active.

    Examples:

        Logout with confirmation prompt:

        >>> import swanlab
        >>> swanlab.logout()

        Force logout without confirmation:

        >>> import swanlab
        >>> swanlab.logout(force=True)
    """
    return raw_logout(force=force)


def raw_logout(force: bool = False) -> bool:
    # 1. 检查是否已登录
    if not apikey.exists():
        console.info("You are not logged in. Use `swanlab login` to login.")
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

    # 3. 重置运行时客户端
    if client.exists():
        client.reset()

    # 4. 删除本地凭证
    apikey.remove()

    console.info("Logout successfully. You can use `swanlab login` to login again.")
    return True
