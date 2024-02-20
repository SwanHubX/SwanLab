#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-20 17:16:44
@File: swanlab\data\modules\video.py
@IDE: vscode
@Description:
    Swanlab的video数据类
"""
from .base import BaseType
from ..utils.file import get_file_hash_numpy_array, get_file_hash_path
import os
from typing import Union, List


class Video(BaseType):
    
    def get_data(self, *args, **kwargs):
        return super().get_data(*args, **kwargs)
    
    
    def get_chart_type(self, *args, **kwargs) -> str:
        return super().get_chart_type(*args, **kwargs)
    
    
    def expect_types(self, *args, **kwargs) -> list:
        return super().expect_types(*args, **kwargs)
    