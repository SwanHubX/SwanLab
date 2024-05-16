#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/16 17:38
@File: use_swanlog.py
@IDE: pycharm
@Description:
    使用swanlog打印日志
"""
from swanlab.log import swanlog

SHOW = "这条日志应该显示"
HIDE = "这条日志不应该显示"

# ---------------------------------- 初始情况 ----------------------------------

swanlog.debug(HIDE)
swanlog.info(SHOW)
swanlog.warning(SHOW)
swanlog.error(SHOW)
swanlog.critical(SHOW)

# ---------------------------------- 安装日志系统 ----------------------------------

swanlog.install(log_level="debug")

swanlog.debug(SHOW)

# ---------------------------------- 卸载日志系统 ----------------------------------

swanlog.uninstall()

swanlog.debug(HIDE)
