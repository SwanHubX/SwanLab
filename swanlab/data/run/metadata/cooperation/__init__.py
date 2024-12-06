"""
@author: cunyue
@file: __init__.py
@time: 2024/12/3 20:17
@description: 合作信息采集
"""

import os
from typing import TypedDict

from swanlab.api import get_http
from swanlab.env import SwanLabEnv, get_mode
from swanlab.package import get_experiment_url
from swanlab.package import get_package_version
from .qing_cloud import get_qing_cloud_info


class SwanLabInfo(TypedDict):
    version: str
    mode: str
    swanlog_dir: str
    exp_url: str
    monitor: int


class CoopInfo(TypedDict):
    swanlab: SwanLabInfo
    qing_cloud: dict


def get_cooperation_info() -> CoopInfo:
    """获取第三方合作信息"""
    qc = get_qing_cloud_info()
    coop = {"swanlab": get_swanlab_info()}
    if qc:
        coop.update({"qing_cloud": qc})
    return coop


def get_swanlab_info() -> SwanLabInfo:
    data: SwanLabInfo = {
        "version": get_package_version(),
        "mode": get_mode(),
        "swanlog_dir": os.environ[SwanLabEnv.SWANLOG_FOLDER.value],
        "monitor": 0,
    }
    try:
        http = get_http()
        data["exp_url"] = get_experiment_url(http.username, http.projname, http.exp_id)
    except ValueError:
        pass
    return data
