"""
@author: caddiesnew
@description: run.fmt 中随 define_scalar 移除而失去调用方的校验函数，v0.11 删除
"""

import warnings

from typing_extensions import deprecated

from swanlab.sdk.internal.pkg import constraints


@deprecated("`safe_validate_chart_name()` has no callers and will be removed in v0.11.")
def safe_validate_chart_name(name):
    """
    检查并清洗图表名称，如果出现非法字符或长度超过限制，返回 None。

    :param name: 待检查的图表名称
    :return: 清洗后的图表名称或 None
    """
    warnings.warn(
        "`safe_validate_chart_name()` has no callers and will be removed in v0.11.",
        FutureWarning,
        stacklevel=2,
    )
    if name is None:
        return None
    try:
        return constraints.ta_chart_name.validate_python(name)
    except Exception:
        return None


@deprecated("`safe_validate_color()` has no callers and will be removed in v0.11.")
def safe_validate_color(color):
    """
    检查并清洗颜色字符串格式，必须是#开头的十六进制颜色代码

    :param color: 待检查的颜色字符串
    :return: 清洗后的颜色字符串或 None
    """
    warnings.warn(
        "`safe_validate_color()` has no callers and will be removed in v0.11.",
        FutureWarning,
        stacklevel=2,
    )
    if color is None:
        return None
    try:
        return constraints.ta_hex_color.validate_python(color)
    except Exception:
        return None
