import atexit
import sys
from ..log import swanlog as swl

# 初始化数据库
from .main import SwanDatabase

swandatabase = SwanDatabase()


# 定义清理函数
def clean_handler():
    if not swl.isError:
        swl.info("train successfully")
        swandatabase.success()
        swl.setSuccess()


# 定义异常处理函数
def except_handler(tp, val, tb):
    swl.error("error happended while training!")
    # 标记实验失败
    swandatabase.fail()
    swl.setError()


# 注册异常处理函数
sys.excepthook = except_handler


# 注册清理函数
atexit.register(clean_handler)

__all__ = ["swandatabase"]
