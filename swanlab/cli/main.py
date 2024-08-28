#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 01:39:04
@File: swanlab\cli\main.py
@IDE: vscode
@Description:
    swanlab脚本命令的主入口
"""
from swanlab.package import get_package_version
import swanlab.cli.commands as C
import click


@click.group(invoke_without_command=True)
@click.version_option(get_package_version(), "--version", "-v", message="SwanLab %(version)s")
def cli():
    pass


# ---------------------------------- 添加其他命令 ----------------------------------


# noinspection PyTypeChecker
cli.add_command(C.login)  # 登录
# noinspection PyTypeChecker
cli.add_command(C.logout)  # 登出

# noinspection PyTypeChecker
cli.add_command(C.watch)  # 启动服务

# noinspection PyTypeChecker
cli.add_command(C.convert)  # 转换命令，用于转换其他实验跟踪工具

# noinspection PyTypeChecker
cli.add_command(C.task)  # 列出、停止、查询任务接口（beta，后续可能删除）

# noinspection PyTypeChecker
cli.add_command(C.launch)  # 启动任务

# noinspection PyTypeChecker
cli.add_command(C.upload)  # 上传文件


if __name__ == "__main__":
    cli()
