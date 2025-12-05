"""
@author: cunyue
@file: log.py
@time: 2025/6/19 16:22
@description: swanlab 终端日志输出器，添加了一些好看的样式
"""

from rich.console import Console
from rich.text import Text

__all__ = ['SwanKitLogger']


class SwanKitLogger:

    def __init__(self, name=__name__.lower(), level: str = "info", file=None):
        self.console = Console(file=file)
        self.__config = {
            "debug": (
                10,
                Text(name, style='grey54 bold', no_wrap=True) + Text(':', style="default"),
            ),
            "info": (
                20,
                Text(name, style='blue bold', no_wrap=True) + Text(':', style="default"),
            ),
            "warning": (
                30,
                Text(name, style='yellow bold', no_wrap=True) + Text(':', style="default"),
            ),
            "error": (
                40,
                Text(name, style='red bold', no_wrap=True) + Text(':', style="default"),
            ),
            "critical": (
                50,
                Text(name, style='red, bold', no_wrap=True) + Text(':', style="default"),
            ),
        }
        self.__can_log = True
        # 默认日志等级为 info
        self.__level: int = self.__config.get(level, self.__config["info"])[0]

    @property
    def level(self):
        """
        获取当前日志等级
        :return:
        """
        for k, v in self.__config.items():
            if v[0] == self.__level:
                return k
        raise AttributeError(f"level {self.__level} not found.")

    @level.setter
    def level(self, level: str):
        """
        设置日志等级
        :param level: 日志等级，可选值为 debug, info, warning, error, critical，如果传入的值不在可选值中，则默认为 info
        """
        if level not in self.__config:
            _level = 20  # info
        else:
            _level = self.__config[level][0]
        self.__level = _level

    def disable_log(self):
        """
        关闭日志输出，实例化时默认开启
        """
        self.__can_log = False

    def enable_log(self):
        """
        开启日志输出
        """
        self.__can_log = True

    def __print(self, log_level: str, *args, **kwargs):
        """
        打印日志
        """
        if not self.__can_log:
            return
        if log_level in self.__config:
            level, prefix = self.__config[log_level]
            if level < self.__level:
                return
            # 处理 sep 参数，即使设置了 sep='' ，前缀后也会有一个空格
            if kwargs.get("sep") == '':
                prefix += " "
            self.console.print(prefix, *args, **kwargs)

    # 发送调试消息
    def debug(self, *args, **kwargs):
        return self.__print("debug", *args, **kwargs)

    # 发送通知
    def info(self, *args, **kwargs):
        return self.__print("info", *args, **kwargs)

    # 发生警告
    def warning(self, *args, **kwargs):
        return self.__print("warning", *args, **kwargs)

    # 发生错误
    def error(self, *args, **kwargs):
        return self.__print("error", *args, **kwargs)

    # 致命错误
    def critical(self, *args, **kwargs):
        return self.__print("critical", *args, **kwargs)
