#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/14 00:51
@File: __init__.py.py
@IDE: pycharm
@Description:
    转换部分，兼容其他可视化工具，转换为swanlab格式
"""


# 使用延迟导入的方式
def __getattr__(name):
    if name == "TFBConverter":
        from .tfb import TFBConverter

        return TFBConverter
    if name == "WandbConverter":
        from .wb import WandbConverter

        return WandbConverter
    
    if name == "MLFlowConverter":
        from .mlf import MLFLowConverter
        
        return MLFLowConverter

    raise AttributeError(f"module 'convert' has no attribute '{name}'")


__all__ = ["TFBConverter", "WandbConverter", "MLFlowConverter"]
