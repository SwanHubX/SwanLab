#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/19 16:46
@File: callback.py
@IDE: pycharm
@Description:
    回调函数注册抽象模块
"""
import atexit
import sys
import traceback

from swanlab.data.run import SwanLabRunState, get_run
from swanlab.log import swanlog
from swanlab.log.type import LogData
from swanlab.swanlab_settings import get_settings
from swanlab.toolkit import SwanKitCallback
from . import utils
from ..store import get_run_store
from ...log.backup import BackupHandler


class SwanLabRunCallback(SwanKitCallback):
    """
    SwanLabRunCallback，回调函数注册类，所有以`on_`和`before_`开头的函数都会在对应的时机被调用
    为了方便管理：
    1. `_`开头的函数为内部函数，不会被调用，且写在最开头
    2. 所有回调按照逻辑上的触发顺序排列
    3. 带有from_*后缀的回调函数代表调用者来自其他地方，比如config、operator等，这将通过settings对象传递
    4. 所有回调不要求全部实现，只需实现需要的回调即可
    """

    def __init__(self):
        self.run_store = get_run_store()
        self.device = BackupHandler()

    def _register_sys_callback(self):
        """
        注册系统回调，内部使用
        """
        sys.excepthook = self._except_handler
        atexit.register(self._clean_handler)

    def _unregister_sys_callback(self):
        """
        注销系统回调，内部使用
        """
        sys.excepthook = sys.__excepthook__
        atexit.unregister(self._clean_handler)

    def _clean_handler(self):
        """
        正常退出清理函数，此函数调用`run.finish`
        """
        run = get_run()
        if run is None:
            return swanlog.debug("SwanLab Runtime has been cleaned manually.")
        # 打印训练结束信息
        utils.print_train_finish(self.run_store.run_name)
        # 如果正在运行
        run.finish() if run.running else swanlog.debug("Duplicate finish, ignore it.")

    @staticmethod
    def _except_handler(tp, val, tb):
        """
        异常退出清理函数
        """
        # 1. 如果是KeyboardInterrupt异常，特殊显示
        if tp == KeyboardInterrupt:
            swanlog.info("KeyboardInterrupt by user")
        else:
            swanlog.info("Error happened while training")
        # 2. 生成错误堆栈
        trace_list = traceback.format_tb(tb)
        error = ""
        for line in trace_list:
            error += line
        error += str(val)
        # 3. 结束运行，注意此时终端错误还没打印
        get_run().finish(SwanLabRunState.CRASHED, error=error)
        assert swanlog.proxied is False, "except_handler should be called after swanlog.stop_proxy()"
        # 4. 打印终端错误，此时终端代理已经停止，不必担心此副作用
        print(error, file=sys.stderr)

    def _terminal_handler(self, log_data: LogData):
        """
        终端输出写入操作
        """
        pass

    def on_run(self, *args, **kwargs):
        # 1. 开启备份
        self.device.start(
            file_dir=self.run_store.file_dir,
            backup_file=self.run_store.backup_file,
            run_name=self.run_store.run_name,
            workspace=self.run_store.workspace,
            visibility=self.run_store.visibility,
            description=self.run_store.description,
            tags=self.run_store.tags,
        )
        # 2. 注册终端输出流代理
        settings = get_settings()
        swanlog.start_proxy(
            proxy_type=settings.log_proxy_type,
            max_log_length=settings.max_log_length,
            handler=self._terminal_handler,
        )
        # 3. 注入系统回调
        self._register_sys_callback()

    def on_stop(self, error: str = None, *args, **kwargs):
        error_epoch = swanlog.epoch + 1
        self.device.stop(error=error, epoch=error_epoch)
        self._unregister_sys_callback()

    def __str__(self):
        raise NotImplementedError("Please implement this method")
