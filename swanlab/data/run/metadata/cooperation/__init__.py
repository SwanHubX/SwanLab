"""
@author: cunyue
@file: __init__.py
@time: 2024/12/3 20:17
@description: 合作信息采集
"""

import os

from swanlab.api import get_http
from swanlab.env import SwanLabEnv
from swanlab.package import get_experiment_url
from swanlab.package import get_package_version
from .qing_cloud import get_qing_cloud_info


def get_cooperation_info():
    qc = get_qing_cloud_info()
    coop = {"swanlab": get_swanlab_info()}
    if qc:
        coop.update({"qing_cloud": qc})
    return coop


def get_swanlab_info():
    data = {
        "version": get_package_version(),
        "mode": os.getenv(SwanLabEnv.MODE.value),
        "swanlog_dir": os.getenv(SwanLabEnv.SWANLOG_FOLDER.value),
    }
    try:
        http = get_http()
        data["exp_url"] = get_experiment_url(http.username, http.projname, http.exp_id)
    except ValueError:
        pass
    return data
