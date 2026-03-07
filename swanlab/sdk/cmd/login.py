"""
@author: cunyue
@file: login.py
@time: 2026/3/6 22:24
@description: swanlab.login 方法，登录到 SwanLab 平台
"""

from typing import Optional


def login(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    web_host: Optional[str] = None,
    save: bool = False,
    timeout: int = 10,
) -> bool:
    """
    Login to SwanLab Cloud.

    This function authenticates your environment with SwanLab. If an API key is
    already configured locally and `relogin` is False, this function will do nothing.

    [Note that] this function should be called before `swanlab.init`.

    :param api_key: str, optional
        Your SwanLab authentication key. If not provided, the SDK will attempt to
        read it from the local environment or prompt you to input it interactively.
    :param relogin: bool, optional
        If True, forces a re-authentication and overwrites the existing API key.
        Defaults to False.
    :param host: str, optional
        The API host URL. If not provided, the default SwanLab cloud host will be used.
    :param web_host: str, optional
        The Web UI host URL, used for printing the dashboard link.
    :param save: bool, optional
        Whether to save the provided api_key to the local netrc/credential file
        for future sessions. Defaults to False.
    :param timeout: int, optional
        Timeout in seconds for the login network request. Defaults to 10.

    :return: bool
        Returns True if login was successful, False otherwise.
    """
    ...
