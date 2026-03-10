"""
@author: cunyue
@file: s3.py
@time: 2026/3/10 14:47
@description: SwanLab S3 文件存储配置
"""

from typing import Optional

from pydantic import BaseModel, Field, field_validator


class S3Settings(BaseModel):
    endpoint: Optional[str] = None

    @field_validator("endpoint", mode="before")
    def validate_endpoint(cls, v: Optional[str]) -> Optional[str]:
        if not v:
            return v

        v = v.strip()
        # 1. 自动剥离用户误填的协议头
        if v.startswith("http://"):
            v = v[7:]
        elif v.startswith("https://"):
            v = v[8:]

        # 2. 剥离末尾可能多余的斜杠
        v = v.rstrip("/")
        return v

    use_ssl: bool = True
    port: Optional[int] = Field(default=None, ge=1, le=65535)
    region: Optional[str] = "us-east-1"
    path_style: bool = False
    bucket: Optional[str] = None
    access_key: Optional[str] = None
    secret_key: Optional[str] = None
