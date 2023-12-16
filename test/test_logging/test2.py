import logging
import logging.config

LOGGING_CONFIG = {
    "version": 1,
    "formatters": {
        "default": {
            "format": "%(asctime)s %(filename)s %(lineno)s %(levelname)s %(message)s",
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
            "filename": "./log.txt",
            "formatter": "default",
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

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger("root")
logger.debug("debug message")
logger.info("info message")
logger.warning("warning message")
logger.error("error message")
logger.critical("critical message")
