#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-20 21:21:45
@File: swanlab\server\controller\experiment.py
@IDE: vscode
@Description:
    实验相关 api 的处理函数
"""

import os
import ujson
import shutil
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409
from fastapi import Request
from urllib.parse import unquote
from ..settings import (
    get_logs_dir,
    get_tag_dir,
    get_files_dir,
    get_exp_dir,
    DB_PATH,
    get_config_path,
)
from ...utils import get_a_lock
from ...utils.file import check_desc_format
import yaml

from ...db import (
    model_talbes,
    connect,
)

from ...db import (
    Project,
    Experiment,
    Tag,
)

__to_list = Experiment.search2list

# ---------------------------------- 通用 ----------------------------------

# 默认项目 id
DEFAULT_PROJECT_ID = Project.DEFAULT_PROJECT_ID
# 实验运行状态
RUNNING_STATUS = Experiment.RUNNING_STATUS


# ---------------------------------- 路由对应的处理函数 ----------------------------------


# 获取实验信息
def get_experiment_info(experiment_id: int):
    """获取实验信息
    1. 数据库中获取实验的基本信息
    2. 从实验目录获取配置信息

    Parameters
    ----------
    experiment_id : int
        实验唯一id
    """

    experiment = Experiment.get(experiment_id).__dict__()
    experiment.pop("project_id")

    # 加载实验配置
    config_path = get_config_path(experiment["run_id"])
    if os.path.exists(config_path):
        with get_a_lock(config_path) as f:
            experiment["config"] = yaml.load(f, Loader=yaml.FullLoader)

    return SUCCESS_200({"experiment": experiment})
