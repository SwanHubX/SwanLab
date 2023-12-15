import logging
from .config import LOGGING_CONFIG
import logging.config


class Swanlog:
    _color_mapping = {
        logging.DEBUG: "\033[36m",  # Cyan
        logging.INFO: "\033[32m",  # Green
        logging.WARNING: "\033[33m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[1;31m",  # Bold Red
    }

    _levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }

    def __init__(self, name=__name__, log_file="app.log", log_level="debug"):
        logging.config.dictConfig(LOGGING_CONFIG)
        self.logger = logging.getLogger(name)
        print(self.logger.handlers)

    def _get_color(self, level):
        # 定义ANSI转义序列
        return self._color_mapping.get(level, "\033[0m")  # Default: Reset color

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
