"""
@author: cunyue
@file: test_utils.py
@time: 2025/9/9 12:32
@description: 测试utils.py中的函数
"""

import os

from swanlab.data.utils import _load_list_from_env


class TestLoadListFromEnv:
    def test_loads_tags_from_env_variable(self):
        os.environ["SWANLAB_TAGS"] = "tag1, tag2, tag3"
        result = _load_list_from_env("SWANLAB_TAGS", None)
        assert result == ["tag1", "tag2", "tag3"]

    def test_returns_value_if_provided(self):
        result = _load_list_from_env("SWANLAB_TAGS", ["custom_tag"])
        assert result == ["custom_tag"]

    def test_returns_none_if_env_variable_not_set(self):
        if "SWANLAB_TAGS" in os.environ:
            del os.environ["SWANLAB_TAGS"]
        result = _load_list_from_env("SWANLAB_TAGS", None)
        assert result is None

    def test_handles_empty_env_variable(self):
        os.environ["SWANLAB_TAGS"] = ""
        result = _load_list_from_env("SWANLAB_TAGS", None)
        assert result == []

    def test_strips_whitespace_from_tags(self):
        os.environ["SWANLAB_TAGS"] = " tag1 , tag2 , tag3 "
        result = _load_list_from_env("SWANLAB_TAGS", None)
        assert result == ["tag1", "tag2", "tag3"]
