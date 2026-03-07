"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:11
@description: SwanLab SDK 版本信息
"""

import json
from pathlib import Path
from typing import Optional

import requests

from swanlab.sdk.pkg.helper import catch_and_return_none

package_path = Path(__file__).resolve().parents[3] / "package.json"


def get_swanlab_version() -> str:
    """获取swanlab的版本号
    :return: swanlab的版本号
    """
    # 读取package.json文件
    with open(package_path, "r") as f:
        return json.load(f)["version"]


@catch_and_return_none()
def get_swanlab_latest_version(timeout=1, url="https://pypi.org/pypi/swanlab/json") -> Optional[str]:
    """
    获取swanlab的最新版本号
    :param timeout: 请求超时时间
    :param url: PyPI API URL
    :return: 最新版本号
    """
    response = requests.get(url, timeout=timeout)
    if response.status_code == 200:
        data = response.json()
        return data["info"]["version"]
    else:
        return None
