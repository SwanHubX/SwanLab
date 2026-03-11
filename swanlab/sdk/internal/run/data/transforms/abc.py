"""
@author: cunyue
@file: abc.py
@time: 2026/3/11 16:32
@description: SwanLab 数据转换模块抽象基类
"""

import inspect
from abc import ABC, abstractmethod
from typing import Any

from google.protobuf.message import Message


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
    >>>     def transform(key: str, step: int, *, data: Any = None, **kwargs: Any) -> Any:
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

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # 只要子类实现了 transform 方法，就检查其签名
        if "transform" in cls.__dict__:
            sig = inspect.signature(cls.transform)
            parameters = list(sig.parameters.values())

            # 至少需要 key 和 step 两个参数
            if len(parameters) < 2:
                raise TypeError(
                    f"The transform method of [{cls.__name__}] must at least include 'key' and 'step' parameters"
                )

            # 1. 检查第一个参数是否为 key
            first_param = parameters[0]
            if first_param.name != "key":
                raise TypeError(f"The first parameter of the transform method of [{cls.__name__}] must be named 'key'")

            # 2. 检查第二个参数是否为 step
            second_param = parameters[1]
            if second_param.name != "step":
                raise TypeError(
                    f"The second parameter of the transform method of [{cls.__name__}] must be named 'step'"
                )

            # 3. 检查从第三个参数开始，是否都是 KEYWORD_ONLY (或者是 **kwargs)
            for param in parameters[2:]:
                # POSITIONAL_OR_KEYWORD 就是普通的参数，如果没有 * 隔离，就会是这种类型
                if param.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
                    raise TypeError(
                        f"The parameters of the transform method of [{cls.__name__}] after 'step' must be keyword-only\n"
                        f"💡 Suggested fix: Add a bare asterisk '*' after the 'step' parameter, for example: def transform(key: str, step: int, *, {param.name}, ...) -> Message:"
                    )

    @staticmethod
    @abstractmethod
    def transform(*args: Any, **kwargs: Any) -> Message:
        """
        将数据转换为Protobuf格式。
        """
        # 套娃加载需求详见 https://github.com/SwanHubX/SwanLab/issues/1367
        ...
