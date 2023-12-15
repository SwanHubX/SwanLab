import logging
from ..env import SWANLAB_CONSOLE_FOLDER
import os
from datetime import datetime

# SWANLAB_CONSOLE_FOLDER 是否存在？如果不存在创建目录
if not os.path.exists(SWANLAB_CONSOLE_FOLDER):
    os.makedirs(SWANLAB_CONSOLE_FOLDER)

LOGGING_CONFIG = {
    "version": 1,
    "formatters": {
        "default": {
            "format": "[%(asctime)s] - [%(filename)s] [%(lineno)s] - [%(name)s] - %(levelname)s %(message)s",
        },
        "detail": {
            "format": "[%(asctime)s]-[%(filename)s]-[%(lineno)s]-[%(name)s] - %(levelname)s: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "DEBUG",  # 输出所有级别的日志到控制台
            "formatter": "default",
        },
        "file": {
            "class": "logging.FileHandler",
            "level": "INFO",  # 输出 INFO 及以上级别的日志到文件
            "filename": os.path.join(SWANLAB_CONSOLE_FOLDER, datetime.now().strftime("%Y-%m-%d") + ".log"),
            "formatter": "detail",
        },
    },
    "loggers": {
        "root": {
            "handlers": ["console", "file"],
            "level": "DEBUG",  # 设置根记录器的级别为 DEBUG，以确保所有级别的日志都会传播到根记录器
        },
    },
    "disable_existing_loggers": False,  # 不禁用现有的记录器
}
