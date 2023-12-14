import sys

log_filename = "output.log"


class Logger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.logfile = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        pass


def init_logger(filename=log_filename):
    # 创建Logger实例，指定保存日志的文件名
    logger = Logger(filename)
    # 将sys.stdout重定向到Logger实例
    sys.stdout = logger
