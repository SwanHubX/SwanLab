#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-05 17:03:54
@File: swanlab/db/table_config.py
@IDE: vscode
@Description:
    数据库表格配置
"""
from .models import (
    Project,
    Experiment,
    Tag,
    Chart,
    Namespace,
    Source,
    Display,
)

tables = [Project, Experiment, Tag, Chart, Namespace, Source, Display]
