from .console import SwanConsoler
from swanlab.utils import FONT
from swankit.log import SwanLabSharedLog


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


class SwanLog(SwanLabSharedLog):

    def __init__(self, name=__name__.lower(), level="info"):
        super().__init__(name=name, level=level)
        # 控制台监控记录器
        self.__consoler = SwanConsoler()
        self.__original_level = level
        self.__installed = False

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
            self.level = log_level
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
        self.level = self.__original_level
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
        return self.__consoler.writer.epoch

    @property
    def file(self):
        if self.__consoler.installed:
            return self.__consoler.writer.file
        else:
            return None

    def can_write(self, level: str) -> bool:
        return self.levels[level] >= self.level

    # 发送调试消息
    @concat_messages
    def debug(self, message):
        super().debug(message)
        return self.can_write("debug")

    # 发送通知
    @concat_messages
    def info(self, message):
        super().info(message)
        return self.can_write("info")

    # 发生警告
    @concat_messages
    def warning(self, message):
        super().warning(message)
        return self.can_write("warning")

    # 发生错误
    @concat_messages
    def error(self, message):
        super().error(message)
        return self.can_write("error")

    # 致命错误
    @concat_messages
    def critical(self, message):
        super().critical(message)
        return self.can_write("critical")

    def reset_console(self):
        """重置控制台记录器"""
        self.__consoler.uninstall()
        # FIXME 这里设置为None，会在test测试中出现ValueError: I/O operation on closed file.
        # self.__consoler = None
