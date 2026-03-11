"""
@author: cunyue
@file: test_scalar.py
@time: 2026/3/11 21:10
@description: 标量处理模块单元测试
"""

import math

import pytest

from swanlab.sdk.internal.run.data.transforms.scalar import Scalar


def test_scalar_instantiation():
    """测试禁止实例化"""
    with pytest.raises(NotImplementedError, match="Scalar should not be instantiated directly"):
        Scalar()


@pytest.mark.parametrize(
    "input_data, expected_number",
    [
        (True, 1.0),  # 布尔值 True
        (False, 0.0),  # 布尔值 False
        (42, 42.0),  # 整型
        (3.1415, 3.1415),  # 浮点型
        ("2.718", 2.718),  # 合法数字字符串
        ("-100", -100.0),  # 负数字符串
    ],
)
def test_scalar_transform_basic_types(input_data, expected_number):
    """测试基础类型及合法字符串的转换"""
    result = Scalar.transform(input_data)
    assert result.number == expected_number


def test_scalar_transform_invalid_string():
    """测试非法字符串转换为 NaN"""
    result = Scalar.transform("invalid_number_string")
    assert math.isnan(result.number)


def test_scalar_transform_numpy():
    """测试真实的 numpy ndarray 转换"""
    import numpy as np

    # 1. 0维标量数组
    res1 = Scalar.transform(np.array(3.14))  # type: ignore
    assert res1.number == pytest.approx(3.14)

    # 2. 1维单元素数组 (numpy 的 .item() 也支持这种)
    res2 = Scalar.transform(np.array([42]))  # type: ignore
    assert res2.number == 42.0

    # 3. 多元素数组拦截
    with pytest.raises(ValueError, match="无法将多元素 Tensor/Array 转换为单一标量记录"):
        Scalar.transform(np.array([1.0, 2.0]))  # type: ignore


def test_scalar_transform_torch():
    """测试真实的 PyTorch Tensor 转换"""
    import torch

    # 1. 0维标量张量
    res1 = Scalar.transform(torch.tensor(2.718))  # type: ignore
    assert res1.number == pytest.approx(2.718)

    # 2. 1维单元素张量
    res2 = Scalar.transform(torch.tensor([-99]))  # type: ignore
    assert res2.number == -99.0

    # 3. 多元素张量拦截
    with pytest.raises(RuntimeError):
        Scalar.transform(torch.tensor([1.0, 2.0, 3.0]))  # type: ignore


def test_scalar_transform_unsupported_type():
    """测试不支持的类型拦截"""
    with pytest.raises(TypeError, match="Unsupported scalar type: list"):
        Scalar.transform([1, 2, 3])  # type: ignore
