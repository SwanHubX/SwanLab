from .utils import init_consoler

# swandatabase对象，使用动态导入的方式有助于环境隔离，比如cli不需要此对象，就不需要导入
_sd = None


def init(experiment_name: str = None, description: str = "", config: dict = {}):
    """初始化swanlab的配置

    Parameters
    ----------
    experiment_name : str, optional
        实验名称, 如果不指定则自动生成一个
    description : str, optional
        实验描述, 如果不指定默认为空
    config : dict, optional
        实验可选配置，在此处可以记录一些实验的超参数等信息
    """
    init_consoler()
    global _sd
    if _sd is not None:
        raise RuntimeError("swanlab has been initialized")
    from .database.main import SwanDatabase

    # 挂载到全局变量
    _sd = SwanDatabase()

    # 注册成功后的清理函数
    import atexit

    atexit.register(_sd.success)

    # 初始化数据库
    _sd.init(
        experiment_name=experiment_name,
        description=description,
        config=config,
    )


def log(data: dict):
    """以字典的形式记录数据，字典的key将作为列名，value将作为记录的值
    例如:
    ```python
    sw.log({"loss": 0.1, "accuracy": 0.9})
    ```
    Parameters
    ----------
    data : dict
        此处填写需要记录的数据
    """
    if _sd is None:
        raise RuntimeError("swanlab has not been initialized")

    if not isinstance(data, dict):
        raise TypeError("log data must be a dict")
    for key in data:
        # 遍历字典的key，记录到本地文件中
        _sd.add(key, data[key])


__all__ = ["init", "log"]
