import sys
import os
from datetime import datetime
from swankit.log import FONT
from swankit.env import create_time


# 检查当前日期是否和控制台日志文件名一致
def check_file_name(func):
    """装饰器，判断是否需要根据日期对控制台输出进行分片存储"""

    def wrapper(self, *args, **kwargs):
        now = datetime.now().strftime("%Y-%m-%d")
        # 检测now是否和self.now一致
        if now != self.now:
            self.now = now
            if hasattr(self, "console") and not self.file.closed:
                self.file.close()
            self.file = open(os.path.join(self.console_folder, self.now + ".log"), "a", encoding="utf-8")
        return func(self, *args, **kwargs)

    return wrapper


class SwanWriterProxy:
    """
    标准输出流拦截代理
    """

    def __init__(self):
        self.epoch = 0
        self.write_callback = None
        self.file = None
        """
        当前正在写入的文件句柄
        """
        self.write_handler = None
        """
        标准输出流原始句柄
        """
        self.__buffer = ""
        """
        上传到云端的缓冲区
        """
        self.console_folder = None
        """
        保存控制台输出文件夹路径
        """
        self.now = None
        """
        当前文件名称（不包含后缀）
        """

    @property
    def can_callback(self) -> bool:
        return self.write_callback is not None

    def set_write_callback(self, func):
        # 封装一层func，加入epoch处理逻辑
        def _func(message):
            self.epoch += 1
            func({"message": message, "create_time": create_time(), "epoch": self.epoch})

        # 封装第二层，加入message处理逻辑以及是否调用逻辑
        def _(message):
            if self.can_callback:
                # 上传到云端
                messages = message.split("\n")
                # 如果长度为1，说明没有换行符
                if len(messages) == 1:
                    self.__buffer = self.__buffer + messages[0]
                # 如果长度大于2，说明其中包含多个换行符
                elif len(messages) > 1:
                    _func(self.__buffer + messages[0])
                    self.__buffer = messages[-1]
                    for m in messages[1:-1]:
                        _func(m)

        self.write_callback = _

    def init(self, path):
        self.console_folder = path
        # path 是否存在
        if not os.path.exists(path):
            os.makedirs(path)
        # 日志文件路径
        self.now = datetime.now().strftime("%Y-%m-%d")
        console_path = os.path.join(path, f"{self.now}.log")
        # 拿到日志文件句柄
        self.file = open(console_path, "a", encoding="utf-8")
        # 封装sys.stdout
        self.write_handler = sys.stdout.write
        a = self.write_handler

        def _(message):
            self.write_handler and self.write_handler(message)
            message = FONT.clear(message)
            self.write_callback and self.write_callback(message)
            self.file.write(message)
            self.file.flush()

        sys.stdout.write = _

    def reset(self):
        sys.stdout.write = self.write_handler
        self.file and self.file.close()
        self.file = None
        self.write_callback = None


class SwanConsoler:
    def __init__(self):
        """
        控制台输出重定向器
        WARNING 一旦此类被初始化，不能将其设置为None，否则会导致输出流无法正常恢复
        """
        self.writer = SwanWriterProxy()
        self.__installed = False

    @property
    def installed(self):
        return self.__installed

    def uninstall(self):
        """重置输出为原本的样子"""
        self.writer.reset()
        self.__installed = False

    def install(self, console_dir):
        """"""
        self.writer.init(console_dir)
        self.__installed = True

    @property
    def write_callback(self):
        return self.writer.write_callback

    def set_write_callback(self, func):
        self.writer.set_write_callback(func)
