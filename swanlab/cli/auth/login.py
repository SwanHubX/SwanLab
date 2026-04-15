"""
@author: caddiesnew
@file: login.py
@time: 2026/4/9 20:38
@description: CLI 登录命令
"""

import click

from swanlab import sdk


@click.command()
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
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host of the swanlab server.",
)
@click.option(
    "--save/--no-save",
    default=True,
    help="Whether to save the API key locally for future sessions.",
)
def login(api_key: str, relogin: bool, host: str, save: bool):
    """Login to the SwanLab cloud."""
    sdk.login_cli(api_key=api_key, relogin=relogin, host=host, save=save)
