"""
@author: cunyue
@file: verify.py
@time: 2025/9/23 11:06
@description: 验证用户当前登录状态
"""

import click
from rich.text import Text

from swanlab.core_python import auth
from swanlab.error import KeyFileError
from swanlab.log import swanlog
from swanlab.package import get_key


@click.command()
def verify():
    """Verify the current login status."""
    try:
        key = get_key()
    except KeyFileError:
        raise click.ClickException("You are not verified. Please use `swanlab login` to login.")
    login_info = auth.code_login(key, save_key=False)
    swanlog.info('You are logged into', login_info.web_host, 'as', Text(login_info.username, "bold"))
