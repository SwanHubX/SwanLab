import sys
import os
from datetime import datetime


class Consoler(sys.stdout.__class__):
    def __init__(self):
        super().__init__(sys.stdout.buffer)

    def init(self, path):
        # 通过当前日期生成日志文件名
        self.now = datetime.now().strftime("%Y-%m-%d")
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
                self.now = now
                if hasattr(self, "console") and not self.console.closed:
                    self.console.close()
                self.console = open(os.path.join(self.console_folder, self.now + ".log"), "a")
            return func(self, *args, **kwargs)

        return wrapper

    @_check_file_name
    def write(self, message):
        self.console.write(message)
        self.console.flush()
        super().write(message)


class SwanConsoler:
    def __init__(self):
        self.consoler = Consoler()

    def init(self, path):
        self.consoler.init(path)
        sys.stdout = self.consoler
