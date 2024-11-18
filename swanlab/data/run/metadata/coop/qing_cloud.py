"""
@author: cunyue
@file: qing_cloud.py
@time: 2024/11/18 15:14
@description: 青云(https://www.qingcloud.com/)元信息采集
"""

import os


BASE_KEYS = ['AICP_PLATFORM', 'AICP_TYPE', 'AICP_NAME', 'AICP_USER_NAME']
RESOURCES_KEYS = [
    'AICP_SPEC_COUNT',
    'AICP_SPEC_GPU',
    'AICP_SPEC_CPU',
    'AICP_SPEC_MEMORY',
    'AICP_SPEC_GPU_NAME',
    'AICP_SPEC_GPU_TYPE',
    'AICP_SPEC_GPU_MEMORY',
    'AICP_HOSTNAME',
    'AICP_HOST_MACHINE',
]


def get_qing_cloud_info():
    plat = os.getenv("AICP_PLATFORM")
    if not plat:
        return None
    return {**get_envs_by_keys(BASE_KEYS), "resources": get_envs_by_keys(RESOURCES_KEYS)}


def get_envs_by_keys(keys: list):
    """
    通过keys获取环境变量，最终返回一个dict，key为keys的值（小写），value为环境变量的值
    """
    return {key.lower(): os.getenv(key) for key in keys}
