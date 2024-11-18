"""
@author: cunyue
@file: official.py
@time: 2024/11/18 15:13
@description: swanlab官方合作信息
"""

import os

from swanlab.api import get_http
from swanlab.data.run.metadata.coop.qing_cloud import get_qing_cloud_info
from swanlab.env import SwanLabEnv
from swanlab.package import get_experiment_url
from swanlab.package import get_package_version


def get_cooperation_info():
    qing_cloud = get_qing_cloud_info()
    coop = {"swanlab": get_swanlab_info()}
    if qing_cloud:
        coop.update({"qing_cloud": qing_cloud})
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
