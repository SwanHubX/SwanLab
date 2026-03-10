"""
@author: cunyue
@file: save.py
@time: 2026/3/10 13:36
@description: SwanLab 文件存储配置
"""

from typing import ClassVar, Type

from pydantic import BaseModel, Field

from swanlab.sdk.internal.settings.s3 import S3Settings


class SaveSettings(BaseModel):
    S3: ClassVar[Type[S3Settings]] = S3Settings

    s3: S3Settings = Field(default_factory=S3Settings)
