"""
@author: cunyue
@file: test_parse.py
@time: 2026/3/14
@description: 测试 config/_parse.py：json_serializable、_adapt_third_party、parse
"""

import argparse
import datetime
from collections import OrderedDict
from dataclasses import dataclass

import pytest

from swanlab.sdk.internal.run.config._parse import _adapt_third_party, json_serializable, parse

# ============================================================
# json_serializable
# ============================================================


class TestJsonSerializable:
    def test_none(self):
        assert json_serializable(None) is None

    def test_bool_true(self):
        result = json_serializable(True)
        assert result is True
        assert type(result) is bool

    def test_bool_false(self):
        result = json_serializable(False)
        assert result is False
        assert type(result) is bool

    def test_int(self):
        assert json_serializable(42) == 42
        assert type(json_serializable(42)) is int

    def test_float(self):
        assert json_serializable(3.14) == 3.14
        assert type(json_serializable(3.14)) is float

    def test_str(self):
        assert json_serializable("hello") == "hello"
        assert type(json_serializable("hello")) is str

    def test_float_nan(self):
        assert json_serializable(float("nan")) == "NaN"

    def test_float_inf(self):
        assert json_serializable(float("inf")) == "Inf"

    def test_float_neg_inf(self):
        assert json_serializable(float("-inf")) == "Inf"

    def test_int_subclass_cast(self):
        """numpy.int64 等子类应转回原生 int"""

        class MyInt(int):
            pass

        result = json_serializable(MyInt(7))
        assert result == 7
        assert type(result) is int

    def test_float_subclass_cast(self):
        class MyFloat(float):
            pass

        result = json_serializable(MyFloat(1.5))
        assert result == 1.5
        assert type(result) is float

    def test_str_subclass_cast(self):
        class MyStr(str):
            pass

        result = json_serializable(MyStr("hi"))
        assert result == "hi"
        assert type(result) is str

    def test_date(self):
        d = datetime.date(2026, 3, 14)
        assert json_serializable(d) == "2026-03-14"

    def test_datetime(self):
        dt = datetime.datetime(2026, 3, 14, 10, 0, 0)
        assert json_serializable(dt) == "2026-03-14T10:00:00"

    def test_list(self):
        assert json_serializable([1, 2.0, "x"]) == [1, 2.0, "x"]

    def test_tuple_to_list(self):
        result = json_serializable((1, 2))
        assert result == [1, 2]
        assert isinstance(result, list)

    def test_nested_list(self):
        assert json_serializable([[1, 2], [3]]) == [[1, 2], [3]]

    def test_dict(self):
        assert json_serializable({"a": 1, "b": "x"}) == {"a": 1, "b": "x"}

    def test_dict_non_str_keys(self):
        result = json_serializable({1: "a"})
        assert result == {"1": "a"}

    def test_mutable_mapping_subclass(self):
        od = OrderedDict([("x", 10), ("y", 20)])
        assert json_serializable(od) == {"x": 10, "y": 20}

    def test_str_fallback(self):
        """无法识别的对象应回退至 str()"""

        class Weird:
            def __str__(self):
                return "weird"

        assert json_serializable(Weird()) == "weird"

    def test_type_error_when_no_str(self):
        """str() 也抛异常时应 raise TypeError"""

        class Bad:
            def __str__(self):
                raise RuntimeError("no str")

        with pytest.raises(TypeError):
            json_serializable(Bad())


# ============================================================
# _adapt_third_party
# ============================================================


class TestAdaptThirdParty:
    def test_argparse_namespace(self):
        ns = argparse.Namespace(lr=0.01, epochs=10)
        result = _adapt_third_party(ns)
        assert result == {"lr": 0.01, "epochs": 10}

    def test_dataclass(self):
        @dataclass
        class Config:
            lr: float
            epochs: int

        result = _adapt_third_party(Config(lr=0.01, epochs=10))
        assert result == {"lr": 0.01, "epochs": 10}

    def test_dataclass_class_itself_raises(self):
        """dataclass 类本身（而非实例）不应被适配"""

        @dataclass
        class Cfg:
            x: int

        with pytest.raises(TypeError):
            _adapt_third_party(Cfg)

    def test_plain_dict_raises(self):
        with pytest.raises(TypeError):
            _adapt_third_party({"a": 1})

    def test_plain_str_raises(self):
        with pytest.raises(TypeError):
            _adapt_third_party("hello")


# ============================================================
# parse
# ============================================================


class TestParse:
    def test_none_returns_empty_dict(self):
        assert parse(None) == {}

    def test_plain_dict(self):
        assert parse({"lr": 0.01, "name": "test"}) == {"lr": 0.01, "name": "test"}

    def test_argparse_namespace(self):
        ns = argparse.Namespace(lr=0.01, epochs=10)
        result = parse(ns)
        assert result == {"lr": 0.01, "epochs": 10}

    def test_dataclass(self):
        @dataclass
        class Cfg:
            lr: float = 0.01
            name: str = "exp"

        result = parse(Cfg())
        assert result == {"lr": 0.01, "name": "exp"}

    def test_dict_with_special_values(self):
        result = parse({"nan": float("nan"), "inf": float("inf"), "ok": 1})
        assert result["nan"] == "NaN"
        assert result["inf"] == "Inf"
        assert result["ok"] == 1

    def test_nested_dict(self):
        result = parse({"nested": {"a": 1}})
        assert result == {"nested": {"a": 1}}

    def test_ordered_dict(self):
        od = OrderedDict([("b", 2), ("a", 1)])
        result = parse(od)
        assert result == {"b": 2, "a": 1}
