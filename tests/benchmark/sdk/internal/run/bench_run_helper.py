"""
@author: cunyue
@file: bench_run_helper.py
@time: 2026/3/12
@description: flatten_dict / validate_key 性能基准测试
"""

from unittest.mock import patch

import pytest

import swanlab.sdk.internal.run.utils_fmt as fmt


@pytest.fixture(params=[10, 100, 1000], ids=["n=10", "n=100", "n=1000"])
def wide_dict(request):
    """参数化宽字典：不同规模的平铺键值对"""
    n = request.param
    return {f"k{i}": i for i in range(n)}


@pytest.fixture
def deep_dict():
    """固定深度 50 的嵌套字典"""
    root: dict = {}
    node = root
    for _ in range(50):
        node["child"] = {}
        node = node["child"]
    node["leaf"] = 1
    return root


# ─────────────────────────── flatten_dict ────────────────────────────


def test_flatten_wide(benchmark, wide_dict):
    """宽字典展开吞吐量（n=10/100/1000）"""
    benchmark(fmt.flatten_dict, wide_dict)


def test_flatten_deep(benchmark, deep_dict):
    """深层嵌套字典展开（depth=50）"""
    benchmark(fmt.flatten_dict, deep_dict)


# ─────────────────────────── validate_key ────────────────────────────


@pytest.fixture(autouse=True)
def reset_warned_keys():
    # noinspection PyProtectedMember
    fmt._WARNED_KEYS.clear()
    yield
    # noinspection PyProtectedMember
    fmt._WARNED_KEYS.clear()


def test_validate_valid(benchmark):
    """合法 key，无需清洗（最快路径）"""
    benchmark(fmt.validate_key, "train/accuracy_top1")


def test_validate_dirty(benchmark):
    """含非法字符的 key，触发 regex 替换"""
    with patch("swanlab.sdk.internal.run.helper.console"):
        benchmark(fmt.validate_key, "  bad key@epoch  ")


def test_validate_long(benchmark):
    """超长 key（300 字符），触发截断"""
    key = "x" * 300
    with patch("swanlab.sdk.internal.run.helper.console"):
        benchmark(fmt.validate_key, key)
