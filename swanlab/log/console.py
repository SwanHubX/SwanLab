import sys
import os
from datetime import datetime
from ..utils import FONT
from swanlab.utils import create_time
from swanlab.utils.judgment import in_jupyter
from io import StringIO


# Consoler 继承的父类
ConsolerParent = sys.stdout.__class__ if not in_jupyter() else StringIO


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


class Consoler(ConsolerParent):
    __init_state = True

    def __init__(self):
        # noinspection PyBroadException
        try:
            # 根据环境进行不同的初始化
            if in_jupyter():
                super().__init__()
            else:
                super().__init__(sys.stdout.buffer)
        except Exception:
            self.__init_state = False
        self.stdout = None
        """
        原始输出流
        """
        self.__buffer = ""
        """
        上传到云端的缓冲区
        """
        self.__epoch = 0
        """
        当前write函数回调调用次数
        """
        self.console_folder = None
        """
        保存控制台输出文件夹路径
        """
        self.now = None
        """
        当前文件名称（不包含后缀）
        """
        self.file = None
        """
        当前正在写入的文件
        """
        # self.__write_callback = None
        # """
        # 注入的上传回调函数本体
        # """
        self.write_callback = None
        """
        封装后的上传回调函数
        """

    @property
    def epoch(self):
        """
        write 函数调用次数（当前行号）
        """
        return self.__epoch

    @property
    def init_state(self) -> bool:
        return self.__init_state

    @property
    def can_callback(self) -> bool:
        return self.write_callback is not None

    def init(self, path, stdout):
        self.console_folder = path
        self.stdout = stdout
        # path 是否存在
        if not os.path.exists(path):
            os.makedirs(path)
        # 日志文件路径
        self.now = datetime.now().strftime("%Y-%m-%d")
        console_path = os.path.join(path, f"{self.now}.log")
        # 如果日志系统初始化失败
        if not self.__init_state:
            with open(console_path, "w", encoding="utf-8") as f:
                f.write("Console Recoder Init Failed!")
        # 拿到日志文件句柄
        self.file = open(console_path, "a", encoding="utf-8")

    def reset(self):
        self.__buffer = ""
        self.__epoch = 0
        self.now = None
        self.file and self.file.close()
        self.file = None
        self.console_folder = None
        self.write_callback = None

    @check_file_name
    def write(self, message):
        """
        重写 write 函数，将输出信息写入到控制台和日志文件中
        """
        # 先写入原始 sys.stdout，即输出到控制台
        self.stdout.write(message)
        self.stdout.flush()
        message = FONT.clear(message)
        self.write_callback and self.write_callback(message)
        self.file.write(message)
        self.file.flush()

    def set_write_callback(self, func):
        # 封装一层func，加入epoch处理逻辑
        def _func(message):
            self.__epoch += 1
            func({"message": message, "create_time": create_time(), "epoch": self.__epoch})

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


class SwanConsoler:
    def __init__(self):
        """
        控制台输出重定向器
        [WARNING] 一旦此类被初始化，不能将其设置为None，否则会导致输出流无法正常恢复
        """
        self.consoler: Consoler = Consoler()
        self.__console_dir = None

    @property
    def installed(self):
        return self.__console_dir is not None

    def uninstall(self):
        """重置输出为原本的样子"""
        if self.installed:
            sys.stdout = sys.__stdout__
            self.consoler.reset()

    def install(self, console_dir):
        self.consoler.init(console_dir, sys.__stdout__)
        self.__console_dir = console_dir
        if self.consoler.init_state:
            sys.stdout = self.consoler

    @property
    def write_callback(self):
        return self.consoler.write_callback

    def set_write_callback(self, func):
        self.consoler.set_write_callback(func)
