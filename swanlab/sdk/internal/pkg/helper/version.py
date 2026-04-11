"""
@author: cunyue
@file: version.py
@time: 2026/4/12 01:19
@description: SwanLab SDK 版本工具
"""

import json
from pathlib import Path
from typing import Optional

import requests

from swanlab.sdk.internal.pkg import safe

package_path = Path(__file__).resolve().parents[4] / "package.json"


def get_swanlab_version() -> str:
    """获取swanlab的版本号
    :return: swanlab的版本号
    """
    # 读取package.json文件
    with open(package_path, "r") as f:
        return json.load(f)["version"]


@safe.decorator(level="debug", message="Failed to fetch swanlab latest version")
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
