"""
@author: cunyue
@file: test_writer.py
@time: 2026/3/14
@description: 测试 config/_writer.py：format/dump 序列化与 write_config 落盘行为
"""

import yaml

from swanlab.sdk.internal.run.components.config.writer import dump_config, format_config, write_config


class TestFormatConfig:
    def test_value_desc_sort_structure(self):
        """每个 key 应归一化为 {value, desc, sort} 结构"""
        content = format_config({"lr": 0.01}, {"lr": 3})

        assert content == {"lr": {"value": 0.01, "desc": "", "sort": 3}}

    def test_sort_from_sort_map(self):
        """sort 字段应取自 sort_map"""
        content = format_config({"a": 1, "b": 2}, {"a": 0, "b": 1})

        assert content["a"]["sort"] == 0
        assert content["b"]["sort"] == 1

    def test_missing_sort_defaults_to_zero(self):
        """sort_map 未包含的 key 默认 sort=0"""
        content = format_config({"x": 99}, {})

        assert content["x"]["sort"] == 0


class TestDumpConfig:
    def test_unicode_allowed(self):
        """中文字符应原样保留，不被转义"""
        text = dump_config(format_config({"名称": "实验一"}, {"名称": 0}))

        assert "名称" in text
        assert "实验一" in text

    def test_round_trip(self):
        """dump 结果可被 yaml.safe_load 完整还原"""
        content = format_config({"lr": 0.01, "epochs": 10, "name": "exp"}, {"lr": 0, "epochs": 1, "name": 2})

        assert yaml.safe_load(dump_config(content)) == content

    def test_empty_config(self):
        """空 config 应序列化为空 YAML 对象（不报错）"""
        text = dump_config(format_config({}, {}))

        assert yaml.safe_load(text) in (None, {})


class TestWriteConfig:
    def test_creates_yaml_file(self, tmp_path):
        path = tmp_path / "config.yaml"
        write_config(path, dump_config(format_config({"lr": 0.01}, {"lr": 0})))

        assert path.exists()

    def test_writes_given_serialized_content(self, tmp_path):
        """写入的应是调用方传入的已序列化文本，不再自行格式化"""
        path = tmp_path / "config.yaml"
        write_config(path, "lr:\n  value: 0.01\n  desc: ''\n  sort: 0\n")

        data = yaml.safe_load(path.read_text())
        assert data == {"lr": {"value": 0.01, "desc": "", "sort": 0}}

    def test_overwrites_existing_file(self, tmp_path):
        """重复调用应全量覆盖，不追加旧内容"""
        path = tmp_path / "config.yaml"
        write_config(path, dump_config(format_config({"lr": 0.01}, {"lr": 0})))
        write_config(path, dump_config(format_config({"epochs": 10}, {"epochs": 0})))

        data = yaml.safe_load(path.read_text())
        assert "lr" not in data
        assert data["epochs"]["value"] == 10
