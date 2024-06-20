def check_library_installed(library_name: str, library_alias: str = None):
    """检查库是否安装，如果未安装则抛出异常"""
    try:
        __import__(library_name)
    except ImportError:
        if library_alias is not None:
            raise RuntimeError(
                f"This contrib module requires {library_alias} to be installed. "
                f"Please install it with command: \n pip install {library_alias}"
            )
        else:
            raise RuntimeError(
                f"This contrib module requires {library_name} to be installed. "
                f"Please install it with command: \n pip install {library_name}"
            )


def check_class_name(class_name):
    """检查类名是否存在"""
    if class_name in globals():
        cls = globals()[class_name]
    else:
        cls = eval(class_name)

    return cls
