"""
@author: cunyue
@file: test_writer.py
@time: 2026/3/14
@description: 测试 config/_writer.py：write_config 序列化与落盘行为
"""

import yaml

from swanlab.sdk.internal.run.config._writer import write_config


class TestWriteConfig:
    def test_creates_yaml_file(self, tmp_path):
        path = tmp_path / "config.yaml"
        write_config(path, {"lr": 0.01}, {"lr": 0})

        assert path.exists()

    def test_value_desc_sort_structure(self, tmp_path):
        """每个 key 应序列化为 {value, desc, sort} 结构"""
        path = tmp_path / "config.yaml"
        write_config(path, {"lr": 0.01}, {"lr": 3})

        data = yaml.safe_load(path.read_text())
        assert data["lr"]["value"] == 0.01
        assert data["lr"]["desc"] == ""
        assert data["lr"]["sort"] == 3

    def test_sort_from_sort_map(self, tmp_path):
        """sort 字段应取自 sort_map"""
        path = tmp_path / "config.yaml"
        write_config(path, {"a": 1, "b": 2}, {"a": 0, "b": 1})

        data = yaml.safe_load(path.read_text())
        assert data["a"]["sort"] == 0
        assert data["b"]["sort"] == 1

    def test_missing_sort_defaults_to_zero(self, tmp_path):
        """sort_map 未包含的 key 默认 sort=0"""
        path = tmp_path / "config.yaml"
        write_config(path, {"x": 99}, {})

        data = yaml.safe_load(path.read_text())
        assert data["x"]["sort"] == 0

    def test_unicode_allowed(self, tmp_path):
        """中文字符应原样保留，不被转义"""
        path = tmp_path / "config.yaml"
        write_config(path, {"名称": "实验一"}, {"名称": 0})

        raw = path.read_text(encoding="utf-8")
        assert "名称" in raw
        assert "实验一" in raw

    def test_empty_config(self, tmp_path):
        """空 config 应写出空 YAML 对象（不报错）"""
        path = tmp_path / "config.yaml"
        write_config(path, {}, {})

        data = yaml.safe_load(path.read_text())
        assert data is None or data == {}

    def test_overwrites_existing_file(self, tmp_path):
        """重复调用应全量覆盖，不追加旧内容"""
        path = tmp_path / "config.yaml"
        write_config(path, {"lr": 0.01}, {"lr": 0})
        write_config(path, {"epochs": 10}, {"epochs": 0})

        data = yaml.safe_load(path.read_text())
        assert "lr" not in data
        assert data["epochs"]["value"] == 10

    def test_multiple_keys(self, tmp_path):
        path = tmp_path / "config.yaml"
        cfg = {"lr": 0.01, "epochs": 10, "name": "exp"}
        sort_map = {"lr": 0, "epochs": 1, "name": 2}
        write_config(path, cfg, sort_map)

        data = yaml.safe_load(path.read_text())
        assert set(data.keys()) == {"lr", "epochs", "name"}
        for k, s in sort_map.items():
            assert data[k]["sort"] == s
