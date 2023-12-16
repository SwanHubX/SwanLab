from ..env import swc
import atexit
import sys

# 注册环境变量，需要在初始化数据库之前注册
swc.init(swc.getcwd(), "train")

# 初始化数据库
from .main import SwanDatabase

swandatabase = SwanDatabase()

# TODO 后续会放到日志对象里
failed = False


# 定义清理函数
def clean_handler():
    print("cleaning...")
    global failed
    if not failed:
        print("success")
        swandatabase.success()


# 定义异常处理函数
def except_handler(tp, val, tb):
    print("error happened")
    global failed
    failed = True
    swandatabase.fail()


# 注册异常处理函数
sys.excepthook = except_handler


# 注册清理函数
atexit.register(clean_handler)

__all__ = ["swandatabase"]
