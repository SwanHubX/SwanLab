from .database import sd


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
    sd.init(
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
    if not isinstance(data, dict):
        raise TypeError("log data must be a dict")
    for key in data:
        # 遍历字典的key，记录到本地文件中
        sd.add(key, data[key])


__all__ = ["init", "log"]
