from .server import SwanWeb
from .database import SwanDatabase
import atexit

swandatabase = SwanDatabase()


init = swandatabase.init

log = swandatabase.add


# 注册退出函数
def exit_handler():
    """退出时执行一些清理操作"""
    swandatabase.success()
    pass


atexit.register(exit_handler)
