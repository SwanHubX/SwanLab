#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:27
@File: __init__.py
@IDE: pycharm
@Description:
    专注于任务启动
"""
import click
import yaml
from .parser import parse
import os

__all__ = ['launch']


@click.command()
@click.option(
    '--file',
    '-f',
    default='swanlab.yml',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help='Designated file to launch',
)
@click.option(
    '--dry-run',
    is_flag=True,
    default=False,
    help='Execute commands without applying changes, only outputting the operations that would be performed.',
)
def launch(file: str, dry_run: bool):
    """
    Launch a task
    """
    file = os.path.abspath(file)
    config = yaml.safe_load(open(file, 'r', encoding='utf-8'))
    if not isinstance(config, dict):
        raise click.FileError(file, hint='Invalid configuration file')
    p = parse(config, file)
    p.parse()
    if dry_run:
        return p.dry_run()
    p.run()
