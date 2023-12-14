import sys
import os
from datetime import datetime
from ..env import SWANLAB_CONSOLE_FOLDER


class Consoler(object):
    def __init__(self, log_directory=SWANLAB_CONSOLE_FOLDER, stream=sys.stdout):
        # 创建保存日志的目录
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)

        # 使用当前日期构建日志文件名
        current_date = datetime.now().strftime("%Y-%m-%d")
        log_filename = os.path.join(log_directory, f"{current_date}.log")

        # 初始化实例
        self.terminal = stream
        self.log = open(log_filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass

    def isatty(self):
        return True


def init_consoler():
    # 创建Logger实例，指定保存日志的文件名
    consoler = Consoler()
    # 将sys.stdout重定向到Logger实例
    sys.stdout = consoler
