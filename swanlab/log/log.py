import logging
from .config import LOGGING_CONFIG
import logging.config
import logging.handlers
from ..env import swc

swc.init(swc.getcwd(), "train")


class Swanlog:
    __color_mapping = {
        logging.DEBUG: "\033[36m",  # Cyan
        logging.INFO: "\033[32m",  # Green
        logging.WARNING: "\033[33m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[1;31m",  # Bold Red
    }

    __levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }

    __status = "running"

    def __init__(self, name=__name__, log_file="output.log", log_level="debug"):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self.__levels[log_level])
        self._create_console_handler()
        self._create_file_handler()

    # 创建控制台记录器
    def _create_console_handler(self, level="debug"):
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        console_handler.setFormatter(formatter)
        console_handler.setLevel(self.__levels[level.lower()])
        self.logger.addHandler(console_handler)

    # 创建日志文件记录器
    def _create_file_handler(self, log_path=None, level="debug"):
        file_handler = logging.FileHandler(swc.output if log_path is None else log_path)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        file_handler.setLevel(self.__levels[level.lower()])
        self.logger.addHandler(file_handler)

    def setOutput(self, log_path=None, level="debug"):
        """
        设置日志文件的存储位置。

        Parameters:
            new_log_path (str): 新的日志文件路径。
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
        if level.lower() in self.__levels:
            console_handler = self.logger.handlers[0]
            console_handler.setLevel(self.__levels[level.lower()])

    def setFileLevel(self, level):
        """
        设置写入日志文件的日志级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        if level.lower() in self.__levels:
            file_handler = self.logger.handlers[1]
            file_handler.setLevel(self.__levels[level.lower()])

    def setLevel(self, level):
        """
        设置日志级别。

        Parameters:
            level (str): 日志级别，可以是 "debug", "info", "warning", "error", 或 "critical".
        """
        if level.lower() in self.__levels:
            self.logger.setLevel(self.__levels[level.lower()])

    def _get_color(self, level):
        # 定义ANSI转义序列
        return self.__color_mapping.get(level, "\033[0m")  # Default: Reset color

    def _reset_color(self):
        # 重置ANSI转义序列
        return "\033[0m"

    def _format_message(self, message, level):
        # 格式化日志消息并添加颜色
        color = self._get_color(level)
        reset_color = self._reset_color()
        return f"{color}{message}{reset_color}"

    def debug(self, message):
        formatted_message = self._format_message(message, logging.DEBUG)
        self.logger.debug(formatted_message)

    def info(self, message):
        formatted_message = self._format_message(message, logging.INFO)
        self.logger.info(formatted_message)

    def warning(self, message):
        formatted_message = self._format_message(message, logging.WARNING)
        self.logger.warning(formatted_message)

    def error(self, message):
        formatted_message = self._format_message(message, logging.ERROR)
        self.logger.error(formatted_message)

    def critical(self, message):
        formatted_message = self._format_message(message, logging.CRITICAL)
        self.logger.critical(formatted_message)
