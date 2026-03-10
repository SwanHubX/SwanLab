"""
@author: cunyue
@file: test_settings_save.py
@time: 2026/3/10 16:30
@description: 测试 SwanLab 文件存储及 S3 配置
"""

import pytest
from pydantic import ValidationError

from swanlab.sdk.internal.settings import Settings


def test_s3_endpoint_validation():
    """
    测试 S3 endpoint 的校验逻辑：
    1. 协议头报错
    2. 末尾斜杠剥离
    3. 空白字符过滤
    """
    # 1. 正常的 endpoint
    s = Settings.model_validate({"save": {"s3": {"endpoint": "custom-s3.example.com"}}})
    assert s.save.s3.endpoint == "custom-s3.example.com"

    # 2. 带有末尾斜杠和多余空格的 endpoint (自动剥离)
    s = Settings.model_validate({"save": {"s3": {"endpoint": "  custom-s3.example.com/  "}}})
    assert s.save.s3.endpoint == "custom-s3.example.com"

    # 3. 空值 (不报错)
    s = Settings.model_validate({"save": {"s3": {"endpoint": ""}}})
    assert s.save.s3.endpoint == ""

    # 4. 包含 http:// 协议头，应当报错
    with pytest.raises(ValueError, match="Invalid s3_endpoint: 'http://custom-s3.example.com'"):
        Settings.model_validate({"save": {"s3": {"endpoint": "http://custom-s3.example.com"}}})

    # 5. 包含 https:// 协议头，应当报错
    with pytest.raises(ValueError, match="Invalid s3_endpoint: 'https://custom-s3.example.com'"):
        Settings.model_validate({"save": {"s3": {"endpoint": "https://custom-s3.example.com"}}})


def test_s3_settings_env_parse(monkeypatch):
    """
    测试通过环境变量解析 S3 相关配置
    根据 SettingsConfigDict，前缀为 SWANLAB_，嵌套分隔符为 _
    """
    # 模拟用户在环境变量中配置私有化 MinIO 存储
    monkeypatch.setenv("SWANLAB_SAVE_S3_ENDPOINT", "s3.my-domain.com")
    monkeypatch.setenv("SWANLAB_SAVE_S3_BUCKET", "my-bucket")
    monkeypatch.setenv("SWANLAB_SAVE_S3_USE_SSL", "false")
    monkeypatch.setenv("SWANLAB_SAVE_S3_PORT", "9000")
    monkeypatch.setenv("SWANLAB_SAVE_S3_PATH_STYLE", "true")
    monkeypatch.setenv("SWANLAB_SAVE_S3_REGION", "cn-north-1")
    monkeypatch.setenv("SWANLAB_SAVE_S3_ACCESS_KEY", "minioadmin")
    monkeypatch.setenv("SWANLAB_SAVE_S3_SECRET_KEY", "minioadmin")

    s = Settings()

    # 验证解析结果
    assert s.save.s3.endpoint == "s3.my-domain.com"
    assert s.save.s3.bucket == "my-bucket"
    assert s.save.s3.use_ssl is False
    assert s.save.s3.port == 9000
    assert s.save.s3.path_style is True
    assert s.save.s3.region == "cn-north-1"
    assert s.save.s3.access_key == "minioadmin"
    assert s.save.s3.secret_key == "minioadmin"


def test_s3_port_boundary():
    """
    测试 S3 port 的边界值验证 (ge=1, le=65535)
    """
    # 1. 正常端口
    s = Settings.model_validate({"save": {"s3": {"port": 443}}})
    assert s.save.s3.port == 443

    # 2. 越界异常 (过小)
    with pytest.raises(ValidationError, match="Input should be greater than or equal to 1"):
        Settings.model_validate({"save": {"s3": {"port": 0}}})

    # 3. 越界异常 (过大)
    with pytest.raises(ValidationError, match="Input should be less than or equal to 65535"):
        Settings.model_validate({"save": {"s3": {"port": 65536}}})


def test_s3_default_values(monkeypatch):
    """
    测试未传入任何参数时，S3 的默认值情况
    需要确保测试环境干净，清理可能干扰的环境变量
    """
    monkeypatch.delenv("SWANLAB_SAVE_S3_ENDPOINT", raising=False)
    monkeypatch.delenv("SWANLAB_SAVE_S3_USE_SSL", raising=False)

    s = Settings()

    assert s.save.s3.endpoint is None
    assert s.save.s3.use_ssl is True
    assert s.save.s3.port is None
    assert s.save.s3.region == "us-east-1"
    assert s.save.s3.path_style is False
    assert s.save.s3.bucket is None
