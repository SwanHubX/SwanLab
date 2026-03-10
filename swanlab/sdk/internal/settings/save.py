"""
@author: cunyue
@file: save.py
@time: 2026/3/10 13:36
@description: SwanLab 文件存储配置
"""

from typing import ClassVar, Type

from pydantic import BaseModel, Field, model_validator

from swanlab.sdk.internal.settings.s3 import S3Settings


class SaveSettings(BaseModel):
    S3: ClassVar[Type[S3Settings]] = S3Settings

    s3: S3Settings = Field(default_factory=S3Settings)

    @model_validator(mode="before")
    @classmethod
    def assemble_nested_env(cls, data: dict) -> dict:
        """
        拦截 Pydantic 因 max_split=1 截断生成的平铺环境变量，
        将其重新组装为嵌套字典，以适配内部结构。
        """
        if isinstance(data, dict):
            # 处理 s3_xxx -> s3: {xxx: ...}
            s3_data = data.get("s3", {})
            if isinstance(s3_data, dict):
                has_update = False
                for key in list(data.keys()):
                    if key.startswith("s3_"):
                        s3_data[key[3:]] = data.pop(key)
                        has_update = True
                if has_update:
                    data["s3"] = s3_data
        return data
