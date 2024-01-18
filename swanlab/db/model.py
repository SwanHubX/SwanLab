#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 13:11:59
@File: swanlab/db/model.py
@IDE: vscode
@Description:
    在此处定义基础模型类
"""
from peewee import Model
from playhouse.shortcuts import model_to_dict
from .settings import swandb


class SwanModel(Model):
    """基础模型类，用于定义数据库表的基本信息"""

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

    def result_to_list(func):
        """将结果转换为列表"""

        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            return [model_to_dict(row) for row in result]

        return wrapper

    class Meta:
        """设置元信息"""

        database = swandb
