"""
@author: cunyue
@file: test_init_load_config.py
@time: 2026/3/10 16:53
@description: 测试 init 中的 load_config 功能
"""

import json
from unittest.mock import MagicMock

import pytest
import yaml

from swanlab.sdk.cmd.init import load_config
from swanlab.sdk.internal.settings import Settings


class TestLoadConfig:
    @pytest.fixture
    def mock_settings(self):
        """模拟带有 run.config 属性的 Settings 对象"""
        settings = MagicMock(spec=Settings)
        settings.run = MagicMock()
        settings.run.config = {"default": "value"}
        return settings

    def test_load_config_from_dict(self, mock_settings):
        """测试：直接传入字典"""
        input_dict = {"key": "value"}
        result = load_config(mock_settings, input_dict)
        assert result == input_dict

    def test_load_config_fallback_to_settings(self, mock_settings):
        """测试：传入 None 时，回退到 run_settings.run.config"""
        result = load_config(mock_settings, None)
        assert result == {"default": "value"}

    def test_load_config_both_none(self):
        """测试：传入 None 且 settings 中也没有配置时"""
        settings = MagicMock()
        settings.run.config = None
        result = load_config(settings, None)
        assert result == {}

    def test_load_config_json_file(self, mock_settings, tmp_path):
        """测试：加载有效的 JSON 文件"""
        config_data = {"api": "v1", "debug": True}
        json_file = tmp_path / "config.json"
        json_file.write_text(json.dumps(config_data))

        result = load_config(mock_settings, str(json_file))
        assert result == config_data

    def test_load_config_yaml_file(self, mock_settings, tmp_path):
        """测试：加载有效的 YAML 文件"""
        config_data = {"project": "swanlab", "tags": ["test", "ml"]}
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(yaml.dump(config_data))

        result = load_config(mock_settings, yaml_file)  # 传入 Path 对象
        assert result == config_data

    def test_load_config_no_extension_json(self, mock_settings, tmp_path):
        """测试：无后缀文件，内容为 JSON"""
        config_data = {"type": "no_ext_json"}
        file = tmp_path / "config_file"
        file.write_text(json.dumps(config_data))

        result = load_config(mock_settings, str(file))
        assert result == config_data

    def test_load_config_no_extension_yaml(self, mock_settings, tmp_path):
        """测试：无后缀文件，内容为 YAML"""
        config_data = {"type": "no_ext_yaml"}
        file = tmp_path / "config_file"
        file.write_text(yaml.dump(config_data))

        result = load_config(mock_settings, str(file))
        assert result == config_data

    def test_load_config_file_not_found(self, mock_settings):
        """测试：文件不存在的情况"""
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            load_config(mock_settings, "non_existent_file.json")

    def test_load_config_invalid_content(self, mock_settings, tmp_path):
        """测试：文件内容格式错误"""
        invalid_file = tmp_path / "bad.json"
        invalid_file.write_text("not a json or yaml {")

        with pytest.raises(ValueError, match="Error parsing config file"):
            load_config(mock_settings, str(invalid_file))

    def test_load_config_invalid_type(self, mock_settings):
        """测试：传入不支持的类型（如整数）"""
        with pytest.raises(ValueError, match="Invalid config type"):
            # noinspection PyTypeChecker
            load_config(mock_settings, 123)  # type: ignore
