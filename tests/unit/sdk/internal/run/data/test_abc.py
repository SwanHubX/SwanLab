"""
@author: cunyue
@file: test_abc.py
@time: 2026/3/11 19:08
@description: SwanLab SDK 数据转换模块抽象基类测试
"""

from typing import Any

import pytest
from google.protobuf.message import Message

# 假设你的代码保存在 swanlab.abc (根据实际路径调整)
from swanlab.sdk.internal.run.data.transforms.abc import TransformType


def test_transform_signature_valid():
    """测试合规的子类签名，不应抛出任何异常"""
    try:

        class ValidTransform1(TransformType):
            @staticmethod
            def transform(key: str, step: int, *, data: Any = None, extra: str = "") -> Message:
                return Message()

        class ValidTransform2(TransformType):
            @staticmethod
            def transform(key: str, step: int, **kwargs: Any) -> Message:
                return Message()
    except Exception as e:
        pytest.fail(f"Compliant subclass definitions should not raise exceptions, but one was raised: {e}")


def test_transform_signature_no_params():
    """测试 transform 没有任何参数（或参数少于2个）的情况"""
    with pytest.raises(TypeError, match=r"must at least include 'key' and 'step' parameters"):

        class InvalidTransformNoParams(TransformType):
            @staticmethod
            def transform(key: str) -> Message:  # 只有一个参数也不行
                return Message()


def test_transform_signature_wrong_first_param():
    """测试 transform 第一个参数名不是 key 的情况"""
    with pytest.raises(TypeError, match=r"first parameter .* must be named 'key'"):

        class InvalidTransformWrongFirst(TransformType):
            @staticmethod
            def transform(wrong_name: str, step: int, *, data: Any) -> Message:
                return Message()


def test_transform_signature_wrong_second_param():
    """测试 transform 第二个参数名不是 step 的情况"""
    with pytest.raises(TypeError, match=r"second parameter .* must be named 'step'"):

        class InvalidTransformWrongSecond(TransformType):
            @staticmethod
            def transform(key: str, wrong_name: int, *, data: Any) -> Message:
                return Message()


def test_transform_signature_missing_asterisk():
    """测试 transform 漏写了 '*' 号，导致后续参数变成位置参数的情况"""
    with pytest.raises(TypeError, match=r"after 'step' must be keyword-only"):

        class InvalidTransformPositional(TransformType):
            @staticmethod
            # 漏掉了 `*` 号，data 变成了普通参数 (POSITIONAL_OR_KEYWORD)
            def transform(key: str, step: int, data: Any) -> Message:
                return Message()


def test_transform_signature_mixed_invalid():
    """测试 transform 中不仅缺少 '*'，还有多个错误位置参数的情况"""
    with pytest.raises(TypeError, match=r"after 'step' must be keyword-only"):

        class InvalidTransformMixed(TransformType):
            @staticmethod
            def transform(key: str, step: int, data: Any, extra: Any) -> Message:
                return Message()


def test_transform_not_implemented_at_creation():
    """
    测试如果子类暂时没有实现 transform 方法，定义时不应报错。
    (因为在抽象类实例化前允许只定义而不实现，报错交由 abc 模块在实例化时接管)
    """
    try:

        class IncompleteTransform(TransformType):
            pass
    except Exception as e:
        pytest.fail(
            f"Abstract subclasses of transform that are not implemented should not raise errors during definition, but threw: {e}"
        )
