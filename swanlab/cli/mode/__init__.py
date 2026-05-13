"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:10
@description: CLI Mode 模块：设置默认运行模式 (disabled / local / online / offline)
"""

import click

from swanlab.sdk import Settings, pkg
from swanlab.sdk.internal.pkg import console


# noinspection PyShadowingNames
def _save_settings(settings: Settings, local: bool):
    """
    将设置保存到文件中
    """
    if local:
        pwd = Settings.get_pwd_config_dir()
        pkg.helper.mkdir_and_append_gitignore(pwd)
        settings.save_to_yaml(pwd, "mode")
    else:
        user = Settings.get_user_config_dir()
        pkg.fs.safe_mkdir(user)
        settings.save_to_yaml(user, "mode")


# noinspection PyShadowingNames
@click.command()
@click.option("--local", is_flag=True, help="Only disable SwanLab for current directory.")
def disabled(local: bool):
    """Disable SwanLab."""
    _save_settings(settings=Settings(mode="disabled"), local=local)

    if local:
        console.info("SwanLab disabled for this directory.")
    else:
        console.info("SwanLab disabled globally.")


# noinspection PyShadowingNames
@click.command()
@click.option("--local", is_flag=True, help="Only set local mode for current directory.")
def local(local: bool):
    """Use local mode for SwanLab."""
    _save_settings(settings=Settings(mode="local"), local=local)

    if local:
        console.info("SwanLab local. Running your script from this directory will save data locally only.")
    else:
        console.info("SwanLab local. Running your script will save data locally only.")


# noinspection PyShadowingNames
@click.command()
@click.option("--local", is_flag=True, help="Only set online mode for current directory.")
def online(local: bool):
    """Use online mode for SwanLab."""
    _save_settings(settings=Settings(mode="online"), local=local)

    if local:
        console.info("SwanLab online. Running your script from this directory will now sync to the cloud.")
    else:
        console.info("SwanLab online. Running your script will now sync to the cloud.")


# noinspection PyShadowingNames
@click.command()
@click.option("--local", is_flag=True, help="Only set offline mode for current directory.")
def offline(local: bool):
    """Use offline mode for SwanLab."""
    _save_settings(settings=Settings(mode="offline"), local=local)

    if local:
        console.info(
            "SwanLab offline. Running your script from this directory will only write metadata locally. "
            "Use `swanlab disabled --local` to completely turn off SwanLab for this directory."
        )
    else:
        console.info(
            "SwanLab offline. Running your script will only write metadata locally. "
            "Use `swanlab disabled` to completely turn off SwanLab."
        )


__all__ = ["disabled", "local", "online", "offline"]
