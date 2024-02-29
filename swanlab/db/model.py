#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-18 13:11:59
@File: swanlab/db/model.py
@IDE: vscode
@Description:
    在此处定义基础模型类
"""
from peewee import Model, OperationalError
from playhouse.shortcuts import model_to_dict
from ..utils import create_time
from .error import NotExistedError
import json


class SwanModel(Model):
    """基础模型类，用于定义数据库表的基本信息"""

    @staticmethod
    def search2dict(result) -> dict:
        """将select、filter的结果转换为字典"""

        dicts = [model_to_dict(row) for row in result]
        if len(dicts) == 1:
            return dicts[0]
        else:
            return dicts

    @staticmethod
    def search2list(result) -> list:
        """将select、filter的结果转换为列表"""

        return [model_to_dict(row) for row in result]

    @staticmethod
    def json2dict(json_str: str) -> dict:
        """将json字符串转换为字典
        如果字符串为空，返回空字典

        Parameters
        ----------
        json_str : str
            传入的json格式字符串，或者空字符串

        Returns
        ----------
        dict :
            转换后的字典

        Raises
        ----------
        TypeError :
            无法转换为字典
        """
        if json_str == "":
            return {}
        try:
            return json.loads(json_str)
        except:
            raise TypeError

    @staticmethod
    def dict2json(dict: dict) -> str:
        """将字典转换为json字符串,如果已经是json字符串，则直接返回

        Parameters
        ----------
        dict : dict
            输入的字典，如果dict为None，则解析为空字符串

        Returns
        -------
        str
            转换后的字符串

                Raises
        ----------

        TypeError :
            无法转换为Json字符串
        """
        if dict is None:
            return ""
        if isinstance(dict, str):
            try:
                if not dict:
                    return ""
                json.loads(dict)
                return dict
            except:
                raise TypeError
        try:
            return json.dumps(dict)
        except:
            raise TypeError

    def save(self, *args, **kwargs):
        """保存数据
        经过swanmodel的覆写，save方法自动更新update_time字段
        """
        self.update_time = create_time()
        super().save(*args, **kwargs)

    @classmethod
    def insert_many(cls, rows, fields=None):
        """批量插入数据
        经过swanmodel的覆写，insert_many方法自动创建create_time和update_time字段
        """
        current_time = create_time()
        rows = [{**row, "create_time": current_time, "update_time": current_time} for row in rows]
        return super().insert_many(rows, fields=fields)

    @classmethod
    def create(cls, **query):
        current_time = create_time()
        query["create_time"] = current_time
        query["update_time"] = current_time
        return super().create(**query)

    @classmethod
    def field_exists(cls, field: str) -> bool:
        """
        判断某个字段是否存在于表中

        Parameters
        ----------
        field : str
            字段名
        """
        # 进行一次查询
        try:
            a = cls.select().where(cls.__getattribute__(cls, field) is not None)
            cls.search2dict(a)
        except OperationalError:
            return False
        return True
