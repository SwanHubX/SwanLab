def init(experiment_name: str = None, description: str = "", config: dict = {}):
    """初始化swanlab的配置
    Parameters
    ----------
    experiment_name : str, optional
        _description_, by default None
    description : str, optional
        _description_, by default ""
    config : dict, optional
        _description_, by default {}
    """
    from .database import sd

    sd.init()
    print("swanlab init")


def log(data: dict):
    """记录数据

    Parameters
    ----------
    data : dict
        _description_
    """
    from .database import sd


__all__ = ["init", "log"]
