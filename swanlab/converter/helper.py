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
