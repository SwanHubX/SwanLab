import re

import pytest

from swanlab.sdk.utils.experiment import (
    generate_color,
    generate_id,
    generate_name,
)

# 定义测试用的精简假数据
MOCK_COLORS = ["#color1", "#color2", "#color3"]
MOCK_ANIMALS = ["mockcat", "mockdog"]
MOCK_ADJECTIVES = ["mockgood", "mocknice"]


class TestGenerateId:
    def test_default_length(self):
        """测试默认生成的 ID 长度"""
        run_id = generate_id()
        assert len(run_id) == 8
        assert run_id.lower() == run_id
        assert run_id.isalnum()
        assert run_id.islower()

    @pytest.mark.parametrize("length", [1, 16, 64])
    def test_custom_length(self, length):
        """测试自定义边界与合法长度"""
        run_id = generate_id(length)
        assert len(run_id) == length

    @pytest.mark.parametrize("invalid_length", [-1, 0, 65, 100])
    def test_invalid_length(self, invalid_length):
        """测试不合法的长度会触发 ValueError"""
        with pytest.raises(ValueError, match="Length must be between 1 and 64."):
            generate_id(invalid_length)


class TestGenerateColor:
    @pytest.fixture(autouse=True)
    def mock_constants(self, monkeypatch):
        """自动在每个测试前替换常量列表"""
        monkeypatch.setattr("swanlab.sdk.utils.experiment.PRESET_COLORS", MOCK_COLORS)

    def test_default_color(self):
        """测试 slug=None 时，返回 Mock 列表中的随机颜色"""
        color = generate_color()
        assert color in MOCK_COLORS

    def test_beauty_color(self):
        """测试 slug='beauty' 时，生成合法的十六进制 HEX 颜色代码"""
        color = generate_color("beauty")
        assert isinstance(color, str)
        assert re.match(r"^#[0-9a-fA-F]{6}$", color), f"Invalid color format: {color}"

    @pytest.mark.parametrize("slug_int", [0, 1, 2, 4, 10])
    def test_int_color(self, slug_int):
        """测试传入 int 时，基于 Mock 列表长度取模返回颜色"""
        expected_color = MOCK_COLORS[slug_int % len(MOCK_COLORS)]
        assert generate_color(slug_int) == expected_color

    def test_fallback_color(self):
        """测试传入不支持的类型或字符串时的兜底逻辑"""
        assert generate_color("invalid_slug") == "#000000"  # type: ignore


class TestGenerateName:
    @pytest.fixture(autouse=True)
    def mock_constants(self, monkeypatch):
        """自动在每个测试前替换常量列表"""
        monkeypatch.setattr("swanlab.sdk.utils.experiment.PRESET_ANIMALS", MOCK_ANIMALS)
        monkeypatch.setattr("swanlab.sdk.utils.experiment.BEAUTY_ADJECTIVES", MOCK_ADJECTIVES)

    def test_default_name(self):
        """测试 slug=None 时，生成 '动物-4位随机字符' 的格式"""
        name = generate_name()
        parts = name.split("-")

        assert len(parts) == 2
        assert parts[0] in MOCK_ANIMALS
        assert len(parts[1]) == 4
        assert parts[1].isalnum()

    def test_beauty_name(self):
        """测试 slug='beauty' 时，生成 '形容词-动物-两位数字' 的格式"""
        name = generate_name("beauty")
        parts = name.split("-")

        assert len(parts) == 3
        assert parts[0] in MOCK_ADJECTIVES
        assert parts[1] in MOCK_ANIMALS
        assert parts[2].isdigit()
        assert 10 <= int(parts[2]) <= 99

    @pytest.mark.parametrize("slug_int", [0, 1, 3, 1024])
    def test_int_name(self, slug_int):
        """测试传入 int 时，基于 Mock 列表长度取模并拼接"""
        expected_animal = MOCK_ANIMALS[slug_int % len(MOCK_ANIMALS)]
        expected_name = f"{expected_animal}-{slug_int}"
        assert generate_name(slug_int) == expected_name

    def test_fallback_name(self):
        """测试兜底机制"""
        name = generate_name("invalid_slug")  # type: ignore
        assert name.startswith("unknown-")
        assert len(name.split("-")[1]) == 4
