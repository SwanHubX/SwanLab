"""
@author: cunyue
@file: constants.py
@time: 2026/3/14
@description: Config 协议常量定义

这些常量是 swanlab.config.v1 协议规范的一部分。
ConfigRecord 中的 path 和 format 字段已从协议中移除，
使用以下约定的常量值以减少序列化后的字节流大小。
"""

# ConfigRecord 的固定路径（相对于 run_dir）
CONFIG_FILE_PATH = "files/config.yaml"

# ConfigRecord 的固定格式
CONFIG_FILE_FORMAT = "yaml"
