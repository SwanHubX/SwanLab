#!/usr/bin/env python# -*- coding: utf-8 -*-
r"""
@DATE: 2025-3-18 21:20
@File: swanlab/swanlab_settings.py
@IDE: cursor
@Description:
    SwanLab全局功能开关，用于管理和控制SwanLab的全局设置
"""
import configparser
import json
import os
from pathlib import Path

from swanlab.env import is_windows, SwanLabMode

try:
    from typing import Annotated, Literal, Optional, Any, Union  # Python 3.9+
except ImportError:
    from typing_extensions import Annotated, Literal, Optional, Any, Union  # Python 3.8

from pydantic import BaseModel, ConfigDict, PositiveInt, StrictBool, DirectoryPath, Field, field_validator


class Settings(BaseModel):
    """
    用户设置类，定义SwanLab支持的所有设置项
    """

    # 使用 Pydantic 的冻结功能代替 @dataclass(frozen=True)
    model_config = ConfigDict(frozen=True)

    # ---------------------------------- 元信息采集部分 ----------------------------------
    # 是否开启元数据采集
    metadata_collect: StrictBool = True
    # 是否采集当前系统环境的硬件信息
    collect_hardware: StrictBool = True
    # 是否采集运行时信息（暂不支持更精细设置，有需求可以添加）
    collect_runtime: StrictBool = True
    # 是否需要主动屏蔽启动命令中的 api key 等隐私信息，默认屏蔽
    security_mask: StrictBool = True
    # ---------------------------------- 其他信息采集 ----------------------------------
    # 是否采集python环境信息
    requirements_collect: StrictBool = True
    # 是否采集conda环境信息
    conda_collect: StrictBool = False
    # ---------------------------------- 硬件监控部分 ----------------------------------
    # 是否开启硬件监控，如果元信息的相关采集被关闭，则此项无效
    hardware_monitor: StrictBool = True
    # 磁盘IO监控的路径
    disk_io_dir: DirectoryPath = Field(
        default_factory=lambda: str(Path(os.environ.get("SystemDrive", "C:")).resolve() if is_windows() else Path("/"))
    )
    hardware_interval: Optional[PositiveInt] = Field(
        default=None,
        ge=5,
        description="Hardware monitoring collection interval, in seconds, minimum value is 5 seconds.",
    )
    # ---------------------------------- 日志上传部分 ----------------------------------
    # 是否开启日志备份功能
    backup: StrictBool = True
    # 日志上传间隔
    upload_interval: PositiveInt = 3
    # 终端日志上传单行最大字符数
    max_log_length: int = Field(ge=500, le=4096, default=1024)
    # 终端日志代理类型，"all"、"stdout"、"stderr"、"none"
    log_proxy_type: Literal["all", "stdout", "stderr", "none"] = "all"

    def filter_changed_fields(self):
        """
        筛选出所有发生变化的设置项
        """
        changed = {}
        # 一些特殊类型例如 disk_io_dir 被设置后无法直接比较（因为是 Path 对象），所以需要特殊处理为字符串
        # 在这里直接转换成JSON然后比较即可
        current_dump = json.loads(self.model_dump_json())
        default_dump = json.loads(Settings().model_dump_json())
        for field_name, default_value in default_dump.items():
            current_value = current_dump[field_name]
            if current_value != default_value:
                changed[field_name] = current_value
        return changed


settings = Settings()


def get_settings():
    """获取当前全局设置"""
    global settings
    return settings


def set_settings(new_settings: Settings) -> Settings:
    global settings
    settings = new_settings
    return settings


def reset_settings():
    global settings
    settings = Settings()


class FolderSettings(BaseModel):
    """
    文件夹级别设置，主要用于读取 INI 文件中的配置，仅在少数场景下会用到，由于这部分配置允许用户自行修改，因此需要额外注意校验逻辑
    """

    mode: str = 'cloud'

    @field_validator('mode')
    def validate_mode(cls, v: Any):
        if v not in SwanLabMode.list():
            # 默认返回 cloud 模式
            return 'cloud'
        return v


def read_folder_settings(folder_path: Union[str, Path], filename: str = 'settings') -> FolderSettings:
    """
    读取指定路径的 INI 配置文件，并返回 FolderSettings 对象
    """
    file_path = Path(folder_path, filename)
    try:
        config = configparser.ConfigParser()
        config.read(file_path)
        kwargs = config['default'] if 'default' in config else {}
        return FolderSettings(**kwargs)
    except configparser.ParsingError:
        return FolderSettings()


def write_folder_settings(folder_path: Union[str, Path], default: dict, filename: str = 'settings'):
    """
    将 FolderSettings 对象写入指定路径的 INI 配置文件
    不保证文件上层目录存在，需要调用方自行保证
    """
    file_path = Path(folder_path, filename)
    config = configparser.ConfigParser()
    config['default'] = default
    with open(file_path, 'w') as configfile:
        config.write(configfile)
