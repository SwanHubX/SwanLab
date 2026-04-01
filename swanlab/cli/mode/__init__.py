"""
@author: cunyue
@file: __init__.py
@time: 2026/4/1 16:10
@description: CLI Mode 模块：设置默认运行模式 (disabled / local / online / offline)
"""

import click


@click.command()
def disabled():
    """Disable SwanLab."""
    # TODO: 接入重构后的 mode 设置逻辑
    pass


@click.command()
def local():
    """Use local mode for SwanLab."""
    # TODO: 接入重构后的 mode 设置逻辑
    pass


@click.command()
def online():
    """Use cloud mode for SwanLab."""
    # TODO: 接入重构后的 mode 设置逻辑
    pass


@click.command()
def offline():
    """Use offline mode for SwanLab."""
    # TODO: 接入重构后的 mode 设置逻辑
    pass


__all__ = ["disabled", "local", "online", "offline"]
