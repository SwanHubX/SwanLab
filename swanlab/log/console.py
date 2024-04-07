import re
import sys
import os
from datetime import datetime
from ..utils import FONT
from swanlab.cloud.files_types import FileType
from typing import Optional


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
        pattern = re.compile(r"swanlab:\s")
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


# 检测是否在 notebook 环境中
def in_notebook():
    try:
        # notebook 中会有 __IPYTHON__，而正常环境没有定义，所以 try
        # 'type: ignore': 可以让 pylance 忽略对变量定义的检查
        __IPYTHON__  # type: ignore
        return True
    except NameError:
        return False


# Consoler 继承的父类
def __consoler_class():
    # 如果在 notebook 中，使用 io.StringIO
    if in_notebook():
        from io import StringIO

        return StringIO
    # 正常环境使用标准输出
    else:
        return sys.stdout.__class__


class Consoler(__consoler_class(), LeverCtl):
    # 记录日志行数
    __sum = 0

    # 上一次输入到标准输出流的消息，用于判断前一个是否换行，即当前行是否为新一行的开头
    # 只有第一次输入时为None
    __previous_message = None

    __init_status = True

    def __init__(self, *args, **kwargs):
        # 根据环境进行不同的初始化
        try:
            if in_notebook():
                super().__init__()
            else:
                super().__init__(sys.stdout.buffer)
        except:
            self.__init_status = False
        self.original_stdout = sys.stdout  # 保存原始的 sys.stdout
        self.pool = None

    @property
    def init_status(self) -> bool:
        return self.__init_status

    def init(self, path):
        # 通过当前日期生成日志文件名
        self.now = datetime.now().strftime("%Y-%m-%d")
        self.console_folder = path
        # path 是否存在
        if not os.path.exists(path):
            os.makedirs(path)
        # 日志文件路径
        console_path = os.path.join(path, f"{self.now}.log")
        # 如果日志系统初始化失败
        if not self.__init_status:
            with open(console_path, "w", encoding="utf-8") as f:
                f.write("Console Recoder Init Failed!")
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
        # 如果直接就是换行
        elif message == "\n":
            # 上一个消息也是以 \n 结尾
            if self.__previous_message.endswith("\n"):
                self.__sum += 1
                message = str(self.__sum) + " \n"
        # 如果在字符串含有 \n
        elif "\n" in message:
            # 通过 \n 切割字符串
            messages = FONT.clear(message).split("\n")
            for index, msg in enumerate(messages):
                # 两种情况不需要处理
                # 为第一个子串，且上一条消息不是以换行结尾，那么不需要添加行号
                # 为最后一个子串，且该子串为空，那么不需要添加行号
                if (index == 0 and not self.__previous_message.endswith("\n")) or (
                    index == len(messages) - 1 and msg == ""
                ):
                    pass
                # 该字串需要单独一行展示，则添加行号
                else:
                    self.__sum += 1
                    msg = str(self.__sum) + " " + msg
                # 如果最后为空，说明原串以 \n 结尾，略过
                if index == len(messages) - 1 and msg == "":
                    pass
                else:
                    self.console.write(msg + "\n")
                    self.upload_message(" ".join(msg.split(' ')[1:]) + "\n")
            self.__previous_message = FONT.clear(message)
            return self.console.flush()
        # 如果是一个头尾不带换行的字符串，需要判断一下前一个message是否带有换行
        else:
            if self.__previous_message.endswith("\n"):
                self.__sum += 1
                message = str(self.__sum) + " " + FONT.clear(message)
            else:
                message = FONT.clear(message)
        self.__previous_message = FONT.clear(message)
        self.console.write(message)
        self.upload_message(" ".join(message.split(' ')[1:]))
        self.console.flush()

    def setLevel(self, level):
        return super().setLevel(level)

    def get_sum(self):
        return self.__sum

    def upload_message(self, message):
        if self.pool is None:
            return
        self.pool.queue.put((FileType.LOG, [message]))


class SwanConsoler:
    def __init__(self):
        self.consoler: Consoler = Consoler()

    def init(self, path):
        self.consoler.init(path)
        if self.consoler.init_status:
            sys.stdout = self.consoler

    def reset(self):
        """重置输出为原本的样子"""
        sys.stdout = self.consoler.original_stdout

    def set_level(self, level):
        """设置控制台打印时，对 swanlog 的收集等级

        Parameters
        ----------
        level : string, optional
            需要设置的等级，["debug", "info", "warning", "error", "critical"] 其中之一
        """
        return self.consoler.setLevel(level)

    def get_level(self):
        return self.consoler.getLevel()

# if __name__ == "__main__":
