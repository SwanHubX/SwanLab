import re
import sys
import os
from datetime import datetime
from ..utils import FONT


class LeverCtl(object):
    # swanlog 的记录权重
    __level_weight = {
        "debug": 10,
        "info": 20,
        "warning": 30,
        "error": 40,
        "critical": 50,
    }

    # 当前记录等级
    __level = "debug"

    def __init__(self, level="debug"):
        """默认收集全部等级的"""

        self.__level = level

    def setLevel(self, level):
        """设置记录等级

        Parameters
        ----------
        level : str, optional
            需要设置的等级，["debug", "info", "warning", "error", "critical"] 其中之一
        """
        self.__level = level

    def __get_weight(self, level):
        """获取记录权重

        Parameters
        ----------
        level : str, optional
            需要获取的等级，["debug", "info", "warning", "error", "critical"] 其中之一

        Returns
        -------
        int
            等级对应权重
        """
        if level in self.__level_weight:
            return self.__level_weight[level]
        else:
            return 0

    def __is_swanlog(self, message):
        """判断信息是否是 swanlog 类的行为

        Parameters
        ----------
        message : string
            打印的信息
        """
        pattern = re.compile(r".*\[SwanLab-(DEBUG|INFO|WARNING|ERROR|CRITICAL)\]:\s+")
        match = pattern.match(message)

        if match:
            # 返回符合要求的部分的小写形式
            return match.group(1).lower()
        else:
            # 不符合要求时返回 None
            return False

    def checkLevel(self, message):
        """检查当前记录等级是否大于等于指定等级

        Parameters
        ----------
        message : string
            打印的信息
        level : str, optional
            需要判断的等级，["debug", "info", "warning", "error", "critical"] 其中之一

        Returns
        -------
        bool
            当前记录等级是否大于等于指定等级
        """
        result = self.__is_swanlog(message)

        if result:
            # 如果匹配到了，进行等级的判断决定是否写入
            return self.__get_weight(result) >= self.__get_weight(self.__level)
        else:
            # 没有匹配到，说明是用户自定义打印，需要记录
            return True

    def getLevel(self):
        """获取当前记录等级

        Returns
        -------
        str
            当前记录等级
        """
        return self.__level


class Consoler(sys.stdout.__class__, LeverCtl):
    # 记录日志行数
    __sum = 0

    # 上一次输入到标准输出流的消息，用于判断前一个是否换行，即当前行是否为新一行的开头
    # 只有第一次输入时为None
    __previous_message = None

    def __init__(self):
        super().__init__(sys.stdout.buffer)
        self.original_stdout = sys.stdout  # 保存原始的 sys.stdout

    def init(self, path, swanlog_level="debug"):
        # 通过当前日期生成日志文件名
        self.now = datetime.now().strftime("%Y-%m-%d")
        self.console_folder = path
        # path 是否存在
        if not os.path.exists(path):
            os.makedirs(path)
        # 日志文件路径
        console_path = os.path.join(path, f"{self.now}.log")
        # 日志文件
        self.console = open(console_path, "a", encoding="utf-8")

    # 检查当前日期是否和控制台日志文件名一致
    def __check_file_name(func):
        """装饰器，判断是否需要根据日期对控制台输出进行分片存储"""

        def wrapper(self, *args, **kwargs):
            now = datetime.now().strftime("%Y-%m-%d")
            # 检测now是否和self.now一致
            if now != self.now:
                self.now = now
                if hasattr(self, "console") and not self.console.closed:
                    self.console.close()
                self.console = open(os.path.join(self.console_folder, self.now + ".log"), "a", encoding="utf-8")
            return func(self, *args, **kwargs)

        return wrapper

    @__check_file_name
    def write(self, message):
        self.original_stdout.write(message)  # 先写入原始 sys.stdout，即输出到控制台
        self.original_stdout.flush()

        # 检查记录等级，高于或等于写入等级即可写入日志文件
        if not self.checkLevel(message):
            return

        if self.__previous_message is None:
            self.__sum += 1
            message = str(self.__sum) + " " + FONT.clear(message)
            self.__previous_message = ""
        # 如果直接就是换行，什么都不做
        elif message == "\n":
            pass
        # 如果以换行符结尾，那么是一个正常的 print 打印，正常处理和输出，单独占一行
        elif message.endswith("\n"):
            if self.__previous_message.endswith("\n"):
                self.__sum += 1
                message = str(self.__sum) + " " + FONT.clear(message)
            else:
                message = FONT.clear(message)
        # 如果是一个不带换行的字符串，需要判断一下前一个message是否带有换行
        else:
            if self.__previous_message.endswith("\n"):
                self.__sum += 1
                message = str(self.__sum) + " " + FONT.clear(message)
            else:
                message = FONT.clear(message)
        self.__previous_message = FONT.clear(message)
        self.console.write(message)
        self.console.flush()

    def setLevel(self, level):
        return super().setLevel(level)

    def getSum(self):
        return self.__sum


class SwanConsoler:
    def __init__(self):
        self.consoler: Consoler = Consoler()

    def init(self, path):
        self.consoler.init(path)
        sys.stdout = self.consoler

    def reset(self):
        """重置输出为原本的样子"""
        sys.stdout = self.consoler.original_stdout

    def setLevel(self, level):
        """设置控制台打印时，对 swanlog 的收集等级

        Parameters
        ----------
        level : string, optional
            需要设置的等级，["debug", "info", "warning", "error", "critical"] 其中之一
        """
        return self.consoler.setLevel(level)

    def getLevel(self):
        return self.consoler.getLevel()


# if __name__ == "__main__":
