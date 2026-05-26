import swanlab


def extract_args(args, kwargs, param_names):
    """从 args/kwargs 按参数名顺序提取值。

    Returns:
        tuple: 按 param_names 顺序返回提取的参数值
    """
    values = []
    for i, name in enumerate(param_names):
        if len(args) > i:
            values.append(args[i])
        else:
            values.append(kwargs.get(name, None))
    return tuple(values)


def ensure_swanlab_init(**init_kwargs):
    """如果 swanlab 尚未初始化，则调用 swanlab.init(**init_kwargs)。"""
    if not swanlab.has_run():
        swanlab.init(**init_kwargs)
