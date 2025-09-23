"""
@author: cunyue
@file: __init__.py
@time: 2025/9/23 17:00
@description: swanlab mode command，等价于自动填写 swanlab.init(mode=xxx)
具体方案是在 swanlog 文件夹下创建 settings 文件，依据 INI 格式存储配置信息，这部分业务存在swanlab.swanlab_settings中
"""

import click

from swanlab.env import create_swanlog_dir
from swanlab.swanlab_settings import write_folder_settings


@click.command()
def disabled():
    """Disable SwanLab"""
    logdir = create_swanlog_dir()
    write_folder_settings(logdir, default={'mode': 'disabled'})


@click.command()
def local():
    """Use local mode for SwanLab"""
    logdir = create_swanlog_dir()
    write_folder_settings(logdir, default={'mode': 'local'})


@click.command()
def online():
    """Use cloud mode for SwanLab"""
    logdir = create_swanlog_dir()
    write_folder_settings(logdir, default={'mode': 'cloud'})


@click.command()
def offline():
    """Use offline mode for SwanLab"""
    logdir = create_swanlog_dir()
    write_folder_settings(logdir, default={'mode': 'offline'})


__all__ = ['disabled', 'local', 'online', 'offline']
