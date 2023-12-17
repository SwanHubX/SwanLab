import atexit, sys, traceback, os
from datetime import datetime
from ..env import swc
from ..log import swanlog as swl

# 注册环境变量，需要在初始化数据库之前注册
swc.init(swc.getcwd(), "train")
# 初始化数据库
from .main import SwanDatabase

swandatabase = SwanDatabase()


# 定义清理函数
def clean_handler():
    if not swl.isError:
        swl.info("train successfully")
        swandatabase.success()
        swl.setSuccess()
        swl.reset_console()


# 定义异常处理函数
def except_handler(tp, val, tb):
    swl.error("Error happended while training, SwanLab will throw it")
    # 标记实验失败
    swandatabase.fail()
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


# 注册异常处理函数
sys.excepthook = except_handler


# 注册清理函数
atexit.register(clean_handler)

__all__ = ["swandatabase"]
