"""
@author: cunyue
@file: abc.py
@time: 2026/3/11 16:32
@description: SwanLab 数据转换模块抽象基类
"""

from abc import ABC, abstractmethod
from typing import Any


class TransformType(ABC):
    """
    SwanLab 数据转换模块抽象基类，定义了将数据转换为Protobuf格式的方法。

    我们约定“TransformType”的相关实现子类实例仅承担“数据传输”职责，数据转换交由每个“TransformType”的实现子类的静态方法“transform”完成。

    设计为静态方法的另一个原因是考虑到大规模数据传输时，静态方法可以避免实例化开销，提高性能。

    每个子类在实现 __init__ 时必须考虑套娃问题，我们约定外层参数的优先级高于内层参数，例如：
    >>> class MyTransform(TransformType):
    >>>     def __init__(self, text: str | MyTransform, foo: Any = None):
    >>>         attrs = self._unwrap(text)
    >>>         self.text = attrs.get("text", text)
    >>>         self.foo = foo if foo is not None else attrs.get("foo")
    >>>     @staticmethod
    >>>     def transform(key:str, data: Any) -> Any:
    >>>         return "example"

    那么此时：

    >>> # 如果内层已经定义了 foo="old"，但外层显式给了 foo="new"
    >>> t = MyTransform(MyTransform("hello", foo="old"), foo="new")
    >>> # 结果应该是 t.foo == "new"
    """

    @classmethod
    def _unwrap(cls, instance_or_val: Any) -> dict:
        """通用的解包辅助函数，提取实例属性，用于实现套娃加载"""
        if isinstance(instance_or_val, TransformType):
            # 返回实例的所有属性（注意排除私有属性）
            return {k: v for k, v in vars(instance_or_val).items() if not k.startswith("_")}
        return {}

    @staticmethod
    @abstractmethod
    def transform(key: str, data: Any) -> Any:
        """
        将数据转换为Protobuf格式，
        此方法必须为静态方法，且不依赖于实例状态，这是为了方便实现“套娃加载”、“懒加载”等策略
        :param key: 数据的键名，用于标识数据类型，此为必填项
        :param data: 待处理的数据
        """
        # 套娃加载需求详见 https://github.com/SwanHubX/SwanLab/issues/1367
        ...
