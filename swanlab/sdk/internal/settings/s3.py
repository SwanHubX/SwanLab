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
        # 1. 发现协议头直接报错，并指导用户如何正确配置
        if v.startswith(("http://", "https://")):
            raise ValueError(
                f"Invalid s3_endpoint: '{v}'. "
                "The endpoint should not contain the protocol (http:// or https://). "
                "Please remove the protocol prefix and use the 's3_use_ssl' configuration instead."
            )

        # 2. 剥离末尾可能多余的斜杠（这个属于无害的格式纠正，可以保留）
        v = v.rstrip("/")
        return v

    use_ssl: bool = True
    port: Optional[int] = Field(default=None, ge=1, le=65535)
    region: Optional[str] = "us-east-1"
    path_style: bool = False
    bucket: Optional[str] = None
    access_key: Optional[str] = None
    secret_key: Optional[str] = None
