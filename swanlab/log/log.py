import logging
import logging.config
import logging.handlers
import sys
from .console import SwanConsoler
from swanlab.utils import FONT


# logging打印格式化类，只负责控制台的相关打印的格式化
class ColoredFormatter(logging.Formatter, FONT):
    def __init__(self, fmt=None, datefmt=None, style="%", handle=None):
        super().__init__(fmt, datefmt, style)
        self.__handle = handle
        # 打印等级对应的颜色装载器
        self.__color_map = {
            logging.DEBUG: self.grey,
            logging.INFO: self.bold_blue,
            logging.WARNING: self.yellow,
            logging.ERROR: self.red,
            logging.CRITICAL: self.bold_red,
        }

    def bold_red(self, s: str) -> str:
        """在终端中加粗的红色字符串

        Parameters
        ----------
        s : str
            需要加粗的字符串

        Returns
        -------
        str
            加粗后的字符串
        """
        # ANSI 转义码用于在终端中改变文本样式
        return self.bold(self.red(s))

    def bold_blue(self, s: str) -> str:
        """在终端中加粗的蓝色字符串

        Parameters
        ----------
        s : str
            需要加粗的字符串

        Returns
        -------
        str
            加粗后的字符串
        """
        return self.bold(self.blue(s))

    def __get_colored_str(self, levelno, message):
        """获取使用打印等级对应的颜色装载的字符串

        Parameters
        ----------
        levelno : logging.levelno
            logging 等级对象
        message : string
            需要装载的颜色
        """

        return self.__color_map.get(levelno)(message)

    def format(self, record):
        """格式化打印字符串
            1. 分割消息头和消息体
            2. 消息头根据 logging 等级装载颜色
            3. 使用空格填充，统一消息头长度为 20 个字符
            4.. 拼接消息头和消息体

        Parameters
        ----------
        record : logging.record
            logging 信息实例

        Returns
        -------
        string
            格式化后的字符串
        """
        log_message = super().format(record)
        self.__handle(log_message + "\n") if self.__handle else None
        # 分割消息，分别处理头尾
        messages: list = log_message.split(":", 1)
        # 填充空格，统一消息头的长度
        message_header = messages[0]
        return f"{self.__get_colored_str(record.levelno, message_header)}:{messages[1]}"


def concat_messages(func):
    """
    装饰器，当传递打印信息有多个时，拼接为一个，并且拦截记录它们
    """

    def wrapper(self, *args, **kwargs):
        # 拼接消息，首先将所有参数转换为字符串，然后拼接
        args = [str(arg) for arg in args]
        message = " ".join(args)
        can_write = func(self, message, **kwargs)
        if can_write and self.file:
            message = self.prefix + FONT.clear(message) + '\n'
            self.file.write(message)
            self.file.flush()
            self.write_callback and self.write_callback(message)

    return wrapper


class SwanLog:
    # 日志系统支持的输出等级
    __levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }

    def __init__(self, name=__name__.lower(), level="info"):
        super().__init__()
        self.prefix = name + ':'
        self.__logger = logging.getLogger(name)
        self.__original_level = self._get_level(level)
        self.__installed = False
        self.__logger.setLevel(self.__original_level)
        # 初始化控制台日志处理器
        self.__handler = logging.StreamHandler(sys.stdout)
        # 添加颜色格式化，并在此处设置格式化后的输出流是否可以被其他处理器处理
        colored_formatter = ColoredFormatter("%(name)s: %(message)s")
        self.__handler.setFormatter(colored_formatter)
        self.enable_log()
        # 控制台监控记录器
        self.__consoler = SwanConsoler()

    def disable_log(self):
        self.__logger.removeHandler(self.__handler)

    def enable_log(self):
        self.__logger.addHandler(self.__handler)

    @property
    def installed(self):
        """
        判断是否已经install
        """
        return self.__installed

    def install(self, console_dir: str = None, log_level: str = None) -> "SwanLog":
        """
        初始化安装日志系统，同一实例在没有执行uninstall的情况下，不可重复安装
        功能是开启标准输出流拦截功能，并设置日志等级
        :param console_dir: 控制台日志文件路径文件夹，如果提供，则会将控制台日志记录到对应文件夹，否则不记录，需要保证文件夹存在
        :param log_level: 日志等级，可以是 "debug", "info", "warning", "error", 或 "critical"，默认为info

        :return: SwanLog实例

        :raises: RuntimeError: 已经安装过日志系统
        :raises: KeyError: 无效的日志级别
        :raises: FileNotFoundError: 控制台日志文件夹不存在
        """
        if self.installed:
            raise RuntimeError("SwanLog has been installed")
        # 设置日志等级
        if log_level is not None:
            self.set_level(log_level)
        # 初始化控制台记录器
        if console_dir:
            self.debug("Init consoler to record console log")
            self.__consoler.install(console_dir)
        self.__installed = True
        return self

    def uninstall(self):
        """
        卸载日志系统，卸载后需要重新安装
        在设计上我们并不希望外界乱用这个函数，所以我们不提供外部调用（不在最外层的__all__中）
        此时将卸载标准输出流拦截功能，并重置日志等级为初始化时的等级
        """
        if not self.installed:
            raise RuntimeError("SwanLog has not been installed")
        self.debug("uninstall swanlog, reset consoler")
        self.__logger.setLevel(self.__original_level)
        self.__consoler.uninstall()
        self.__installed = False

    @property
    def write_callback(self):
        return self.__consoler.write_callback

    def set_write_callback(self, func):
        self.__consoler.set_write_callback(func)

    @property
    def epoch(self):
        """
        获取当前日志的 epoch
        """
        return self.__consoler.consoler.epoch

    def set_level(self, level):
        """
        Set the logging level of the logger.

        :param level: The level to set the logger to. This should be one of the following:
            - "debug"
            - "info"
            - "warning"
            - "error"
            - "critical"

        :raises: KeyError: If an invalid level is passed.
        """
        self.__logger.setLevel(self._get_level(level))

    # 获取对应等级的logging对象
    def _get_level(self, level):
        """私有属性，获取等级对应的 logging 对象

        Parameters
        ----------
        level : string
            日志级别，可以是 "debug", "info", "warning", "error", 或 "critical"

        Returns
        -------
        logging.level : object
            logging 模块中的日志等级

        Raises
        ------
        KeyError
            无效的日志级别
        """
        if level.lower() in self.__levels:
            return self.__levels.get(level.lower())
        else:
            raise KeyError("log_level must be one of ['debug', 'info', 'warning', 'error', 'critical']: %s" % level)

    @property
    def file(self):
        if self.__consoler.installed:
            return self.__consoler.consoler.file
        else:
            return None

    def can_write(self, level: str) -> bool:
        return self._get_level(level) >= self.__logger.level

    # 发送调试消息
    @concat_messages
    def debug(self, message):
        self.__logger.debug(message)
        return self.can_write("debug")

    # 发送通知
    @concat_messages
    def info(self, message):
        self.__logger.info(message)
        return self.can_write("info")

    # 发生警告
    @concat_messages
    def warning(self, message):
        self.__logger.warning(message)
        return self.can_write("warning")

    # 发生错误
    @concat_messages
    def error(self, message):
        self.__logger.error(message)
        return self.can_write("error")

    # 致命错误
    @concat_messages
    def critical(self, message):
        self.__logger.critical(message)
        return self.can_write("critical")

    def reset_console(self):
        """重置控制台记录器"""
        self.__consoler.uninstall()
        # FIXME 这里设置为None，会在test测试中出现ValueError: I/O operation on closed file.
        # self.__consoler = None
