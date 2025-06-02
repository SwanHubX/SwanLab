"""
@author: cunyue
@file: local.py
@time: 2025/2/13 14:59
@description: 本地模式下的回调
local模式（目前）将自动调用swanboard，如果不存在则报错
"""

from swankit.callback import ColumnInfo

from swanlab.log.type import LogData
from swanlab.swanlab_settings import get_settings

try:
    # noinspection PyPackageRequirements
    import swanboard
except ImportError:
    raise ImportError("Please install swanboard to use 'local' mode: pip install 'swanlab[dashboard]'")

from importlib.metadata import version

package_version = version("swanboard")
if package_version != "0.1.8b1":
    raise ImportError(
        "Your swanboard version does not match, please use this command to install the matching version: pip install 'swanlab[dashboard]'"
    )


import json
import os
import sys
import traceback
from datetime import datetime
from typing import Union, Tuple, Optional, TextIO

from swankit.callback.models import RuntimeInfo, MetricInfo
from swankit.core import SwanLabSharedSettings

from swanlab.data.run.callback import SwanLabRunCallback
from swanlab.data.run.main import get_run, SwanLabRunState
from swanlab.env import SwanLabEnv
from swanlab.log import swanlog


class LocalRunCallback(SwanLabRunCallback):

    def __init__(self):
        super(LocalRunCallback, self).__init__()
        self.board = swanboard.SwanBoardCallback()
        # 当前日志写入文件的句柄
        self.file: Optional[TextIO] = None

    @staticmethod
    def _traceback_error(tb, val):
        """
        获取traceback信息
        """
        trace_list = traceback.format_tb(tb)
        html = ""
        for line in trace_list:
            html += line
        html += str(val)
        return html

    @staticmethod
    def _error_print(tp):
        """
        错误打印
        """
        # 如果是KeyboardInterrupt异常
        if tp == KeyboardInterrupt:
            swanlog.info("KeyboardInterrupt by user")
        else:
            swanlog.info("Error happened while training")

    @staticmethod
    def _init_logdir(logdir: Union[str, None] = None) -> None:
        """
        根据传入的logdir,初始化日志文件夹
        ---
        Args:
            logdir: 日志文件夹路径

        Return:
            None

        Step:
            1: 参数检查
            2: 环境变量设置
            3: 默认路径
            4: .gitignore
        """
        env_key = SwanLabEnv.SWANLOG_FOLDER.value
        # 如果传入了logdir，则将logdir设置为环境变量，代表日志文件存放的路径
        # 如果没有传入logdir，则使用默认的logdir, 即当前工作目录下的swanlog文件夹，但是需要保证目录存在
        if logdir is None:
            logdir = os.environ.get(env_key) or os.path.join(os.getcwd(), "swanlog")

        logdir = os.path.abspath(logdir)
        try:
            os.makedirs(logdir, exist_ok=True)
            if not os.access(logdir, os.W_OK):
                raise IOError(f"no write permission for path: {logdir}")
        except Exception as error:
            raise IOError(f"Failed to create or access logdir: {logdir}, error: {error}")

        os.environ[env_key] = logdir

        # 如果logdir是空的，创建.gitignore文件，写入*
        if not os.listdir(logdir):
            with open(os.path.join(logdir, ".gitignore"), "w", encoding="utf-8") as f:
                f.write("*")

    def __str__(self):
        return "SwanLabLocalRunCallback"

    def _except_handler(self, tp, val, tb):
        """
        异常处理
        """
        self._error_print(tp)
        # 结束运行
        get_run().finish(SwanLabRunState.CRASHED, error=self._traceback_error(tb, tp(val)))
        if tp != KeyboardInterrupt:
            print(self._traceback_error(tb, tp(val)), file=sys.stderr)

    def _clean_handler(self):
        run = get_run()
        if run is None:
            return swanlog.debug("SwanLab Runtime has been cleaned manually.")
        self._train_finish_print()
        # 如果正在运行
        run.finish() if run.running else swanlog.debug("Duplicate finish, ignore it.")

    def _write_handler(self, log_data: LogData):
        log_name = f"{datetime.now().strftime('%Y-%m-%d')}.log"
        if self.file is None:
            # 如果句柄不存在，则创建
            self.file = open(os.path.join(self.settings.console_dir, log_name), "a", encoding="utf-8")
        elif os.path.basename(self.file.name) != log_name:
            # 如果句柄存在，但是文件名不一样，则关闭句柄，重新打开
            self.file.close()
            self.file = open(os.path.join(self.settings.console_dir, log_name), "a", encoding="utf-8")
        # 写入日志
        for content in log_data["contents"]:
            self.file.write(content['message'] + '\n')
            self.file.flush()

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs):
        self._init_logdir(logdir)
        self.board.on_init(proj_name)

    def before_run(self, settings: SwanLabSharedSettings, *args, **kwargs):
        self.settings = settings
        # path 是否存在
        if not os.path.exists(self.settings.run_dir):
            os.mkdir(self.settings.run_dir)

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        #  FIXME num 在 dashboard 中被要求传递但是没用上 🤡
        self.board.before_init_experiment(run_id, exp_name, description, colors=colors, num=1)

    def on_run(self):
        settings = get_settings()
        if settings.log_proxy_type != "none":
            swanlog.start_proxy(settings.log_proxy_type, -1, self._write_handler)
        # 注入系统回调
        self._register_sys_callback()
        # 打印信息
        self._train_begin_print()
        self._watch_tip_print()

    def on_runtime_info_update(self, r: RuntimeInfo, *args, **kwargs):
        # 更新运行时信息
        if r.requirements is not None:
            r.requirements.write(self.settings.files_dir)
        if r.metadata is not None:
            r.metadata.write(self.settings.files_dir)
        if r.config is not None:
            r.config.write(self.settings.files_dir)

    def on_log(self, *args, **kwargs):
        self.board.on_log(*args, **kwargs)

    def on_column_create(self, column_info: ColumnInfo, *args, **kwargs):
        self.board.on_column_create(column_info)

    def on_metric_create(self, metric_info: MetricInfo, *args, **kwargs):
        # 出现任何错误直接返回
        if metric_info.error:
            return
        # 屏蔽非自定义指标，因为现在本地不支持系统指标
        if metric_info.column_info.cls != "CUSTOM":
            return
        # ---------------------------------- 保存指标数据 ----------------------------------
        self.settings.mkdir(os.path.dirname(metric_info.metric_file_path))
        self.settings.mkdir(os.path.dirname(metric_info.summary_file_path))
        with open(metric_info.summary_file_path, "w+", encoding="utf-8") as f:
            f.write(json.dumps(metric_info.metric_summary, ensure_ascii=False))
        with open(metric_info.metric_file_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(metric_info.metric, ensure_ascii=False) + "\n")

        # ---------------------------------- 保存媒体字节流数据 ----------------------------------
        if metric_info.metric_buffers is None:
            return
        for i, r in enumerate(metric_info.metric_buffers):
            if r is None:
                continue
            # 组合路径
            path = os.path.join(metric_info.swanlab_media_dir, metric_info.column_info.kid)
            os.makedirs(path, exist_ok=True)
            # 写入数据
            with open(os.path.join(path, metric_info.metric["data"][i]), "wb") as f:
                f.write(r.getvalue())

    def on_stop(self, error: str = None, *args, **kwargs):
        """
        训练结束，取消系统回调
        此函数被`run.finish`调用
        """
        # 写入错误信息
        if error is not None:
            with open(self.settings.error_path, "a") as fError:
                print(datetime.now(), file=fError)
                print(error, file=fError)
        # 打印信息
        self._watch_tip_print()
        # 取消注册系统回调
        self._unregister_sys_callback()

        self.board.on_stop()
