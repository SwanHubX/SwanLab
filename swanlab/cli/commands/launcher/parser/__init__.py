#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/8/22 14:29
@File: __init__.py.py
@IDE: pycharm
@Description:
    解析配置文件
"""
from . import v1
from .model import LaunchParser
from typing import Dict, List
import click

__all__ = ['parse']

parsers: Dict[str, List[LaunchParser.__class__]] = {
    'swanlab/v1': v1.parsers
}


def parse(config: dict, path: str) -> LaunchParser:
    version = config.get("apiVersion")
    if not parsers.get(version):
        raise click.UsageError(f"Unknown api version: {version}")
    kind = config.get("kind")
    for parser in parsers[version]:
        if parser.__type__() == kind:
            return parser(config, path)
    raise click.UsageError(f"Unknown kind: {kind}")
