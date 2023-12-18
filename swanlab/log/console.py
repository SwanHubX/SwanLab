import sys
import os
from datetime import datetime


class Consoler(sys.stdout.__class__):
    def __init__(self):
        super().__init__(sys.stdout.buffer)
        self.original_stdout = sys.stdout  # 保存原始的 sys.stdout

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
        self.console = open(console_path, "a", encoding="utf-8")

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
                self.console = open(os.path.join(self.console_folder, self.now + ".log"), "a", encoding="utf-8")
            return func(self, *args, **kwargs)

        return wrapper

    @_check_file_name
    def write(self, message):
        self.console.write(message)
        self.console.flush()
        self.original_stdout.write(message)  # 同时写入原始 sys.stdout
        self.original_stdout.flush()

    @_check_file_name
    def add(self, message: str):
        """此接口用于写入额外的信息到日志文件中，但是不会写入到控制台

        Parameters
        ----------
        message : str
            写入的信息
        """
        self.console.write(message)
        self.console.flush()


class SwanConsoler:
    def __init__(self):
        self.consoler: Consoler = Consoler()
        self.add: function = self.consoler.add

    def init(self, path):
        self.consoler.init(path)
        sys.stdout = self.consoler

    def reset(self):
        """重置输出为原本的样子"""
        sys.stdout = self.consoler.original_stdout
