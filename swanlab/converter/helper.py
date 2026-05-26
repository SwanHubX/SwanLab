import inspect


def extract_args(fn, args, kwargs, param_names):
    """从 args/kwargs 按函数签名提取指定参数值。

    Returns:
        tuple: 按 param_names 顺序返回提取的参数值
    """
    try:
        bound = inspect.signature(fn).bind_partial(*args, **kwargs)
    except (TypeError, ValueError):
        return tuple(kwargs.get(name, None) for name in param_names)

    return tuple(bound.arguments.get(name, None) for name in param_names)
