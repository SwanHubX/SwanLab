"""
@author: cunyue
@file: _parse.py
@time: 2026/3/14
@description: Config 值序列化——纯函数，无副作用，无外部依赖
"""

import argparse
import datetime
import json
import math
from collections.abc import MutableMapping
from dataclasses import asdict, is_dataclass
from typing import Any, cast

__all__ = ["parse"]

# 基础原生类型（子类需要向上转型）
_BASE_TYPES = (int, float, str)


def json_serializable(obj: Any) -> Any:
    """
    将任意 Python 对象递归转换为 JSON 可序列化的基础类型。

    特殊值处理：
    - float NaN  → 字符串 "NaN"
    - float Inf  → 字符串 "Inf"
    - bool 优先于 int 处理（bool 继承自 int）
    - 子类类型（如 numpy.int64）→ 强制转回父类原生类型

    :raises TypeError: 对象无法被序列化
    """
    if obj is None:
        return None

    # float 特殊值须在 _BASE_TYPES 循环前处理
    if type(obj) is float:
        if math.isnan(obj):
            return "NaN"
        if math.isinf(obj):
            return "Inf"

    # bool 继承自 int，必须先判断
    if type(obj) is bool:
        return obj

    for t in _BASE_TYPES:
        if type(obj) is t:
            return obj
        # 子类（如 numpy.float64）→ 转为原生类型
        if isinstance(obj, t):
            return t(obj)

    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()

    if isinstance(obj, (list, tuple)):
        return [json_serializable(item) for item in obj]

    if isinstance(obj, (dict, MutableMapping)):
        return {str(k): json_serializable(v) for k, v in obj.items()}

    try:
        return str(obj)
    except Exception:
        raise TypeError(f"Object {obj!r} is not JSON serializable")


def _adapt_third_party(data: Any) -> dict:
    """
    适配第三方配置对象，转换为普通 dict。

    支持：
    - omegaconf.DictConfig
    - mmengine.Config
    - argparse.Namespace
    - dataclass 实例

    :raises TypeError: 未能命中任何适配器
    """
    # omegaconf
    try:
        import omegaconf  # noqa

        if isinstance(data, omegaconf.DictConfig):
            return omegaconf.OmegaConf.to_container(data, resolve=True, throw_on_missing=True)  # type: ignore
    except ImportError:
        pass

    # mmengine
    try:
        import mmengine  # noqa

        if isinstance(data, mmengine.Config):
            return mmengine.Config.to_dict(data)
    except ImportError:
        pass

    # argparse.Namespace（标准库，无需 try/except）
    if isinstance(data, argparse.Namespace):
        return vars(data)

    # dataclass 实例（注意排除 dataclass 类本身）
    if is_dataclass(data) and not isinstance(data, type):
        # noqa: 虽然警告 'dataclasses.asdict' method should be called on dataclass instances，但此处 data 已经是 dataclass 实例了
        return asdict(cast(Any, data))

    raise TypeError


def parse(config: Any) -> dict:
    """
    将 config 转换为 JSON 可序列化的 dict。

    转换策略（按优先级依次尝试）：
    1. 第三方类型适配（omegaconf、mmengine、argparse、dataclass）
    2. json_serializable 递归转换
    3. json.dumps / json.loads 兜底

    :raises TypeError: 所有策略均失败
    """
    if config is None:
        return {}

    # 1. 第三方类型适配
    try:
        return _adapt_third_party(config)
    except TypeError:
        pass

    # 2. json_serializable
    try:
        result = json_serializable(config)
        if isinstance(result, dict):
            return result
    except TypeError:
        pass

    # 3. JSON round-trip 兜底
    try:
        return json.loads(json.dumps(config))
    except Exception as e:
        raise TypeError(f"config {config!r} is not a JSON-serializable dict: {e}")
