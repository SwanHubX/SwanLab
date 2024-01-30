#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:43:21
@File: swanlab/server/settings.py
@IDE: vscode
@Description:
    服务部分基础配置
"""
import os
import mimetypes
from ..env import get_swanlog_dir

"""
在此处注册静态文件路径，因为静态文件由vue框架编译后生成，在配置中，编译后的文件存储在/swanlab/template中
入口文件为index.html，网页图标为logo.ico，其他文件为assets文件夹中的文件
"""
# 注册静态文件路径
mimetypes.add_type("application/javascript", ".js")
mimetypes.add_type("text/css", ".css")
# 静态文件路径
FILEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(FILEPATH, "template")
ASSETS = os.path.join(TEMPLATE_PATH, "assets")
INDEX = os.path.join(TEMPLATE_PATH, "index.html")
# swanlog文件夹路径
SWANLOG_DIR = get_swanlog_dir()


# ---------------------------------- 实验目录下的各个路径 ----------------------------------

# tags 存储目录
LOGS = "logs"
# files 实验环境目录
FILES = "files"
# console 实验环境目录
CONSOLE = "console"
# media 实验生成媒体目录
MEDIA = "media"
# 实验环境(元数据)
META_DATA = "swanlab-metadata.json"
# 实验依赖
REQUIREMENTS = "requirements.txt"
# 数据库名
DB_NAME = "runs.swanlab"
# 数据库路径
DB_PATH = os.path.join(SWANLOG_DIR, DB_NAME)
# 实验配置
EXPERIMENT_CONFIG = "config.yaml"


def get_exp_dir(name) -> str:
    """实验目录路径"""
    return os.path.join(SWANLOG_DIR, name)


def get_logs_dir(name) -> str:
    """logs 目录路径，存放 tag 数据"""
    return os.path.join(SWANLOG_DIR, name, LOGS)


def get_files_dir(name) -> str:
    """files 目录路径"""
    return os.path.join(SWANLOG_DIR, name, FILES)


def get_console_dir(name) -> str:
    """console 目录路径"""
    return os.path.join(SWANLOG_DIR, name, CONSOLE)


def get_config_path(name) -> str:
    """实验配置存储路径"""
    return os.path.join(get_files_dir(name), EXPERIMENT_CONFIG)


def get_meta_path(name) -> str:
    """实验环境存储路径"""
    return os.path.join(SWANLOG_DIR, name, FILES, META_DATA)


def get_requirements_path(name) -> str:
    """实验依赖存储路径"""
    return os.path.join(SWANLOG_DIR, name, FILES, REQUIREMENTS)


def get_tag_dir(name, tag) -> str:
    """获取 tag 对应的目录路径"""
    return os.path.join(SWANLOG_DIR, name, LOGS, tag)


def get_media_dir(name, tag) -> str:
    """获取 media 对应的目录路径"""
    return os.path.join(SWANLOG_DIR, name, MEDIA, tag)
