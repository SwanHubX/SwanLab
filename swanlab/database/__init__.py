import atexit, sys, traceback, os
from datetime import datetime
from ..env import swc
from ..log import swanlog as swl
from ..utils.file import check_key_format

from .modules import BaseType


sd = None


# 定义清理函数
def clean_handler():
    if not swl.isError:
        swl.info("train successfully")
        sd.success()
        swl.setSuccess()
        swl.reset_console()


# 定义异常处理函数
def except_handler(tp, val, tb):
    swl.error("Error happended while training, SwanLab will throw it")
    # 标记实验失败
    sd.fail()
    swl.setError()
    # 记录异常信息
    # 追踪信息
    traceList = traceback.format_tb(tb)
    html = repr(tp) + "\n"
    html += repr(val) + "\n"
    for line in traceList:
        html += line + "\n"

    if os.path.exists(swc.error):
        swl.warning("Error log file already exists, append error log to it")
    # 写入日志文件
    with open(swc.error, "a") as fError:
        print(datetime.now(), file=fError)
        print(html, file=fError)
    # 重置控制台记录器
    swl.reset_console()
    raise tp(val)


def init(experiment_name: str = None, description: str = "", config: dict = {}, *args, **kwargs):
    """初始化swanlab的配置

    Parameters
    ----------
    experiment_name : str, optional
        实验名称, 如果不指定则自动生成一个
    description : str, optional
        实验描述, 如果不指定默认为空
    config : dict, optional
        实验可选配置，在此处可以记录一些实验的超参数等信息
    kwargs : dict, optional
        其他额外的非必要设置项目
    """
    global sd
    # 注册环境变量，需要在初始化数据库之前注册
    swc.init(swc.getcwd(), "train")
    # 初始化数据库
    from .main import SwanDatabase

    sd = SwanDatabase()

    # 注册异常处理函数
    sys.excepthook = except_handler

    # 注册清理函数
    atexit.register(clean_handler)

    # 初始化数据库
    sd.init(
        experiment_name=experiment_name,
        description=description,
        config=config,
    )
    # 初始化日志对象

    swl.init(swc.output, level=kwargs.get("log_level", "info"))
    swl.debug("SwanLab database initialized")
    swl.debug("Swanlab will take over all the print information of the terminal from now on")
    swl.info("Run data will be saved locally in " + swc.exp_folder)
    swl.info("Experiment_name: " + sd.experiment.name)
    swl.info("Run `swanlab watch` to view SwanLab Experiment Dashboard")


def log(data: dict, step: int = None):
    """以字典的形式记录数据，字典的key将作为列名，value将作为记录的值
    例如:
    ```python
    sw.log({"loss": 0.1, "accuracy": 0.9})
    ```
    Parameters
    ----------
    data : dict
        此处填写需要记录的数据
    step: int
        当前记录的步数，如果不传则默认当前步数为'已添加数据数量+1'
    """
    if sd is None:
        raise RuntimeError("swanlab is not initialized")
    if not isinstance(data, dict):
        raise TypeError("log data must be a dict")
    if step is not None and (not isinstance(step, int) or step < 0):
        raise TypeError("'step' must be an integer not less than zero.")
    for key in data:
        # 遍历字典的key，记录到本地文件中
        d = data[key]
        # 检查key的类型
        check_key_format(key)
        # 检查数据类型，data[key]必须是int，float或者可以被float化的类型，或者swanlab.BaseType的子类
        if not isinstance(data[key], (int, float, BaseType)):
            try:
                d = float(data[key])
            except:
                raise TypeError("log data must be int, float, swanlab.BaseType or can be converted to float")
        # 添加数据
        sd.add(key, d, step=step)


__all__ = ["log", "init", "BaseType"]
