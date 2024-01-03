import logging
import logging.config
import logging.handlers
import sys
from .console import SwanConsoler
from ..utils import FONT


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


# 控制台打印格式化类，只负责控制台的相关打印的格式化
class ColoredFormatter(logging.Formatter, FONT):
    def __init__(self, fmt=None, datefmt=None, style="%", handle=None):
        super().__init__(fmt, datefmt, style)
        self.__handle = handle
        # 打印等级对应的颜色装载器
        self.__color_map = {
            logging.DEBUG: self.grey,
            logging.INFO: self.green,
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
        target_length = 20
        # 填充空格，统一消息头的长度
        message_header = messages[0] + ":" + " " * max(0, target_length - len(messages[0]))
        return f"{self.__get_colored_str(record.levelno, message_header)} {messages[1]}"


# 日志文件格式化类
class CustomFileHandler(logging.FileHandler):
    def emit(self, record):
        # 在写入日志之前把颜色剔除
        processed_message = FONT.clear(record.getMessage())
        record.msg = processed_message
        # 调用父类的 emit 方法写入日志
        super().emit(record)


class SwanLog(Logsys):
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

    def init(
        self,
        path,
        level=None,
        console_level=None,
        file_level=None,
        console_path=None,
    ):
        """初始化内部打印器
            初始化的顺序最好别变，下面的一些设置方法没有使用查找式获取处理器，而是直接用索引获取的
            所以 handlers 列表中，第一个是控制台处理器，第二个是日志文件处理器

        Parameters
        ----------
        path : string
            日志文件的路径，打印会记录到文件
        level : string, optional
            全局日志级别，设置该等级会同时影响控制台和文件, by default None
        console_level : string, optional
            控制台日志级别，仅影响控制台
        file_level : string, optional
            文件日志级别，高于或等于该级别即记录到文件
        console_path: str, optional
            控制台日志文件路径，如果提供，则会将控制台日志记录到文件,否则不记录
        """

        # 初始化控制台记录器
        if self.__consoler is None and console_path:
            self.debug("init consoler")
            self.__consoler = SwanConsoler()
            self.__consoler.init(console_path)

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
        console_handler = logging.StreamHandler(sys.stdout)
        # 添加颜色格式化，并在此处设置格式化后的输出流是否可以被其他处理器处理
        colored_formatter = ColoredFormatter("[%(name)s-%(levelname)s]: %(message)s")
        console_handler.setFormatter(colored_formatter)
        console_handler.setLevel(self._getLevel(level))
        self.logger.addHandler(console_handler)

    # 创建日志文件记录器
    @_check_init
    def _create_file_handler(self, log_path, level="debug"):
        file_handler = CustomFileHandler(log_path, encoding="utf-8")
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

    def setCollectionLevel(self, level):
        """
        设置日志收集级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        self.__consoler.setLevel(level)

    def getCollectionLevel(self):
        """
        获取日志收集级别。

        Returns:
            str: 日志收集级别。
        """
        return self.__consoler.getLevel()

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
