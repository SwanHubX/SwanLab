"""
@author: cunyue
@file: __init__.py
@time: 2026/4/1 16:09
@description: CLI 认证模块：login / logout / verify
"""

import click

from swanlab.sdk.cmd.login import interactive_login


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
    interactive_login(api_key=api_key, relogin=relogin, host=host, save=save)


@click.command()
def logout():
    """Logout from the SwanLab cloud."""
    # TODO: 接入重构后的 logout 逻辑
    pass


@click.command()
def verify():
    """Verify the current login status."""
    # TODO: 接入重构后的 verify 逻辑
    pass
