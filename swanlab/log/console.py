import sys
import os
from datetime import datetime


class SwanConsoler(sys.stdout.__class__):
    def __init__(self):
        super().__init__(sys.stdout.buffer)

    def init(self, path):
        # 通过当前日期生成日志文件名
        # self.now = datetime.now().strftime("%Y-%m-%d")
        self.now = "2023-12-16"
        self.console_folder = path
        # path 是否存在
        if not os.path.exists(path):
            os.makedirs(path)
        # 日志文件路径
        console_path = os.path.join(path, f"{self.now}.log")
        # 日志文件
        self.console = open(console_path, "a")

    # 检查当前日期是否和控制台日志文件名一致
    def _check_file_name(func):
        """装饰器，判断是否需要根据日期对控制台输出进行分片存储"""

        def wrapper(self, *args, **kwargs):
            now = datetime.now().strftime("%Y-%m-%d")
            # 检测now是否和self.now一致
            if now != self.now:
                print("console recoder path changed")
                self.now = now
                self.console = open(os.path.join(self.console_folder, self.now), "a")
            return func(*args, **kwargs)

        return wrapper

    @_check_file_name
    def write(self, message):
        self.console.write(message)
        self.console.flush()
        super().write(message)


def init_consoler(path):
    # 创建Logger实例，指定保存日志的文件名
    consoler = SwanConsoler()
    consoler.init(path)
    # 将sys.stdout重定向到Logger实例
    sys.stdout = consoler
