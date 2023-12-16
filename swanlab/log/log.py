import logging
import logging.config
import logging.handlers


class Logsys:
    # 日志系统状态：running / success / error
    __status = "running"

    def __init__(self):
        self.__status = "running"

    def setSuccess(self):
        if self.isRunning:
            self.__status = "success"
        else:
            raise KeyError("%s is not running" % self.__status)

    def setError(self):
        if self.isRunning:
            self.__status = "error"
        else:
            raise KeyError("%s is not running" % self.__status)

    @property
    def isSuccess(self) -> bool:
        return self.__status == "success"

    @property
    def isError(self) -> bool:
        return self.__status == "error"

    @property
    def isRunning(self) -> bool:
        return self.__status == "running"


class Swanlog(Logsys):
    # 转义之后的颜色系统
    __color_mapping = {
        logging.DEBUG: "\033[36m",  # Cyan
        logging.INFO: "\033[32m",  # Green
        logging.WARNING: "\033[33m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[1;31m",  # Bold Red
    }

    # 日志系统支持的输出等级
    __levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }

    def __init__(self, name=__name__, log_file="output.log", level="debug"):
        super()
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self._getLevel(level))

    def init(self, path):
        # 初始化的顺序最好别变，下面的一些设置方法没有使用查找式获取处理器，而是直接用索引获取的
        # 所以 handlers 列表中，第一个是控制台处理器，第二个是日志文件处理器
        self._create_console_handler()
        self._create_file_handler(path)

    # 创建控制台记录器
    def _create_console_handler(self, level="debug"):
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        console_handler.setFormatter(formatter)
        console_handler.setLevel(self.__levels[level.lower()])
        self.logger.addHandler(console_handler)

    # 创建日志文件记录器
    def _create_file_handler(self, log_path, level="debug"):
        file_handler = logging.FileHandler(log_path)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        file_handler.setLevel(self.__levels[level.lower()])
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

    # 获取对应等级的logging对象
    def _getLevel(self, level):
        if level.lower() in self.__levels:
            return self.__levels.get(level.lower())
        else:
            raise KeyError("Invalid log level: %s" % level)

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

    # 发送调试消息
    def debug(self, message):
        formatted_message = self._format_message(message, logging.DEBUG)
        self.logger.debug(formatted_message)

    # 发送通知
    def info(self, message):
        formatted_message = self._format_message(message, logging.INFO)
        self.logger.info(formatted_message)

    # 发生警告
    def warning(self, message):
        formatted_message = self._format_message(message, logging.WARNING)
        self.logger.warning(formatted_message)

    # 发生错误
    def error(self, message):
        formatted_message = self._format_message(message, logging.ERROR)
        self.logger.error(formatted_message)

    # 致命错误
    def critical(self, message):
        formatted_message = self._format_message(message, logging.CRITICAL)
        self.logger.critical(formatted_message)
