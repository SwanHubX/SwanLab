"""
@author: cunyue
@file: test_mime.py
@time: 2026/7/29
@description: 测试 swanlab.sdk.internal.core_python.pkg.mime 模块

测试分组：
  - TestGuessTypeRegistered : 已注册扩展（.yaml/.yml）的类型推断
  - TestGuessTypePassthrough : 标准扩展的透传
  - TestGuessTypeFallback    : 未知扩展/无扩展的兜底
  - TestGuessTypeInput       : str 与 PathLike 输入兼容
"""

import pytest

from swanlab.sdk.internal.core_python.pkg.mime import guess_type


class TestGuessTypeRegistered:
    @pytest.mark.parametrize("name", ["config.yaml", "a/b/c.yml", "swanlab.yaml", "UPPER.YAML", "Mixed.YmL"])
    def test_yaml_extensions(self, name):
        assert guess_type(name) == "application/yaml"


class TestGuessTypePassthrough:
    def test_known_extension(self):
        # .txt 在所有平台均映射为 text/plain
        assert guess_type("readme.txt") == "text/plain"


class TestGuessTypeFallback:
    @pytest.mark.parametrize("name", ["data.unknownext", "archive.zzz", "noextension", ""])
    def test_unknown_or_missing_extension(self, name):
        assert guess_type(name) == "application/octet-stream"


class TestGuessTypeInput:
    def test_accepts_path(self, tmp_path):
        f = tmp_path / "config.yaml"
        f.write_text("key: value")
        assert guess_type(f) == "application/yaml"

    def test_accepts_str(self):
        assert guess_type("config.yaml") == "application/yaml"
