#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-16 11:00:55
@File: swanlab\db\setting.py
@IDE: vscode
@Description:
    一些配置
"""
import os
from ..env import get_swanlog_dir
from peewee import SqliteDatabase, Model
from playhouse.shortcuts import model_to_dict

# 连接到SQLite数据库，设置数据库文件名为runs.swanlab
swandb = SqliteDatabase(os.path.join(get_swanlog_dir(), "runs.swanlab"))


class SwanModel(Model):
    """基础模型类，用于定义数据库表的基本信息"""

    @staticmethod
    def result_to_dict(func):
        """将结果转换为字典"""

        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            dicts = [model_to_dict(row) for row in result]
            if len(dicts) == 1:
                return dicts[0]
            else:
                return dicts

        return wrapper

    class Meta:
        """设置元信息"""

        database = swandb
