import logging
import logging.config
import logging.handlers
from .console import SwanConsoler
from ..env import swc


class Logsys:
    # 日志系统状态：running / success / error
    __status = "running"

    def __init__(self):
        self.__status = "running"

    def setSuccess(self):
        if self.isRunning:
            self.__status = "success"
        else:
            raise Exception("current status is %s. You can only set success while runnging" % self.__status)

    def setError(self):
        if self.isRunning:
            self.__status = "error"
        else:
            raise Exception("current status is %s. You can only set success while runnging" % self.__status)

    @property
    def isSuccess(self) -> bool:
        return self.__status == "success"

    @property
    def isError(self) -> bool:
        return self.__status == "error"

    @property
    def isRunning(self) -> bool:
        return self.__status == "running"


# 新增的带颜色的格式化类
class ColoredFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style="%", handle=None):
        super().__init__(fmt, datefmt, style)
        self.__handle = handle

    _color_mapping = {
        logging.DEBUG: "\033[90m",  # Grey
        logging.INFO: "\033[32m",  # Green
        logging.WARNING: "\033[33m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[1;31m",  # Bold Red
    }

    def format(self, record):
        log_message = super().format(record)
        self.__handle(log_message + "\n") if self.__handle else None
        color = self._color_mapping.get(record.levelno, "\033[0m")  # Default: Reset color
        reset_color = "\033[0m"
        # 分割消息，分别处理头尾
        messages: list = log_message.split(":", 1)
        target_length = 20
        message_header = messages[0] + ":" + " " * max(0, target_length - len(messages[0]))
        return f"{color}{message_header}{reset_color} {messages[1]}"


class Swanlog(Logsys):
    # 日志系统支持的输出等级
    __levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }

    def __init__(self, name=__name__, level="debug"):
        super()
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self._getLevel(level))
        self.__consoler: SwanConsoler = None

    def init(self, path, level=None, console_level=None, file_level=None):
        # 初始化的顺序最好别变，下面的一些设置方法没有使用查找式获取处理器，而是直接用索引获取的
        # 所以 handlers 列表中，第一个是控制台处理器，第二个是日志文件处理器

        # 初始化控制台记录器
        if self.__consoler is None and swc.isTrain:
            self.debug("init consoler")
            self.__consoler = SwanConsoler()
            self.__consoler.init(swc.console_folder)

        self._create_console_handler()
        self._create_file_handler(path)
        if level:
            self.logger.setLevel(self._getLevel(level))
        if console_level:
            self.setConsoleLevel(console_level)
        if file_level:
            self.setFileLevel(file_level)

    # 检测日志处理器是否重复注册
    def _check_init(func):
        """装饰器，防止多次注册处理器"""

        def wrapper(self, *args, **kwargs):
            if len(self.logger.handlers) == 2:
                return self.debug("init more than once")
            result = func(self, *args, **kwargs)
            return result

        return wrapper

    # 创建控制台记录器
    @_check_init
    def _create_console_handler(self, level="debug"):
        console_handler = logging.StreamHandler()
        handle = None if self.__consoler is None else self.__consoler.add
        # 添加颜色格式化，并在此处设置格式化后的输出流是否可以被其他处理器处理
        colored_formatter = ColoredFormatter("[%(name)s-%(levelname)s]: %(message)s", handle=handle)
        console_handler.setFormatter(colored_formatter)
        console_handler.setLevel(self._getLevel(level))
        self.logger.addHandler(console_handler)

    # 创建日志文件记录器
    @_check_init
    def _create_file_handler(self, log_path, level="debug"):
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        formatter = logging.Formatter("%(name)s %(levelname)s [%(asctime)s] %(message)s")
        file_handler.setFormatter(formatter)
        file_handler.setLevel(self._getLevel(level))
        self.logger.addHandler(file_handler)

    def setOutput(self, log_path=None, level="debug"):
        """
        设置日志文件的存储位置。

        Parameters:
            log_path (str): 日志文件路径。
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        file_handler = self.logger.handlers[1]  # Assuming file handler is the second handler
        self.logger.removeHandler(file_handler)
        self._create_file_handler(log_path, level)

    def setConsoleLevel(self, level):
        """
        设置控制台输出的日志级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        console_handler = self.logger.handlers[0]
        console_handler.setLevel(self._getLevel(level))

    def setFileLevel(self, level):
        """
        设置写入日志文件的日志级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        file_handler = self.logger.handlers[1]
        file_handler.setLevel(self._getLevel(level))

    def setLevel(self, level):
        """
        设置日志级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        self.logger.setLevel(self._getLevel(level))

    # 获取对应等级的logging对象
    def _getLevel(self, level):
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
            raise KeyError(
                "Invalid log level: %s, level must be one of ['debug', 'info', 'warning', 'error', 'critical']" % level
            )

    # 发送调试消息
    def debug(self, message):
        self.logger.debug(message)

    # 发送通知
    def info(self, message):
        self.logger.info(message)

    # 发生警告
    def warning(self, message):
        self.logger.warning(message)

    # 发生错误
    def error(self, message):
        self.logger.error(message)

    # 致命错误
    def critical(self, message):
        self.logger.critical(message)

    def reset_console(self):
        """重置控制台记录器"""
        self.__consoler.reset()
        self.__consoler = None
