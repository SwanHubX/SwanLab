"""
@author: cunyue
@file: __init__.py
@time: 2026/3/12
@description: SwanLab SDK 运行模块，涉及：
1. 数据处理 (基于 Event-Bus 事件驱动架构)
2. 运行、实验上下文管理
3. 触发异步微批处理落盘与回调
"""

import atexit
import signal
import sys
import threading
import traceback
from functools import cached_property, wraps
from pathlib import Path
from types import TracebackType
from typing import Any, Literal, Mapping, Optional, Type, Union, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus import RunEmitter
from swanlab.sdk.internal.bus.events import MetricLogEvent, RunFinishEvent, RunStartEvent, ScalarDefineEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core import create_core
from swanlab.sdk.internal.pkg import console, log
from swanlab.sdk.internal.run.config import (
    Config as _ConfigClass,
)
from swanlab.sdk.internal.run.config import (
    create_run_config,
    create_unbound_run_config,
    deactivate_run_config,
)
from swanlab.sdk.internal.run.transforms import Audio, Image, Text, Video, normalize_media_input
from swanlab.sdk.typings.run import FinishType
from swanlab.sdk.typings.run.column import ScalarXAxisType
from swanlab.sdk.typings.run.transforms import CaptionsType
from swanlab.sdk.typings.run.transforms.audio import AudioDatasType, AudioRatesType
from swanlab.sdk.typings.run.transforms.image import ImageDatasType, ImageFilesType, ImageModesType, ImageSizesType
from swanlab.sdk.typings.run.transforms.text import TextDatasType
from swanlab.sdk.typings.run.transforms.video import VideoDatasType

from . import utils_fmt as fmt
from .consumer import BackgroundConsumer
from .record_builder import RecordBuilder

__all__ = ["Run", "has_run", "get_run", "set_run", "clear_run"]


def with_lock(func):
    """线程安全装饰器，自动获取和释放外部API锁
    Python 3.9 中无法定义在类内部
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        with self._api_lock:
            return func(self, *args, **kwargs)

    return wrapper


def with_run(cmd: str):
    """
    run api装饰器，确保当前 run 实例激活
    """

    def decorator(f):
        @wraps(f)
        def wrapper(self, *args, **kwargs):
            if self._state != "running":
                raise RuntimeError(f"`{cmd}` requires an active Run, call `swanlab.init()` first.")
            return f(self, *args, **kwargs)

        return wrapper

    return decorator


class Run:
    """A SwanLab run for tracking experiments.

    This class represents a single experiment run and provides methods for logging
    metrics, artifacts, and metadata. Typically created via `swanlab.init()` and
    accessed via `swanlab.get_run()`.

    Attributes:

        config: Configuration object for storing hyperparameters.

        id: Unique identifier for this run.

        run_dir: Local directory where run data is stored.

        url: Cloud URL for this run (None in local/offline mode).

        project_url: Cloud URL for the project (None in local/offline mode).

    Examples:

        Access run properties:

        >>> import swanlab
        >>> run = swanlab.init(mode="local", project="my_project")
        >>> print(run.id)
        >>> print(run.run_dir)
        >>> swanlab.finish()

        Log metrics:

        >>> import swanlab
        >>> run = swanlab.init(mode="local")
        >>> run.log({"loss": 0.5, "accuracy": 0.95})
        >>> swanlab.finish()
    """

    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        self._state: Union[FinishType, Literal["running"]] = "running"
        # 外部API锁，防止并发调用
        self._api_lock = threading.RLock()

        # 事件发射器：唯一的队列写入入口，可注入给内部系统组件
        self._emitter = RunEmitter(maxsize=100_000)

        # Core：Record 落盘与后端交互的统一入口
        self._core = create_core(ctx)

        # 记录构建器：负责将事件转换为 Record
        self._builder = RecordBuilder(ctx)

        # 后台消费者：从 emitter.queue 消费事件并落盘
        self._consumer = BackgroundConsumer(ctx, self._emitter.queue, self._builder, self._core)
        self._consumer.start()

        # 触发启动事件
        ts = Timestamp()
        ts.GetCurrentTime()
        self._emitter.emit(RunStartEvent(timestamp=ts))

        # TODO: 硬件监控与metadata采集

        # 绑定 config 模块到运行上下文
        if self._ctx.config.settings.mode != "disabled":
            self.config: _ConfigClass = create_run_config(self._ctx.config_file, self._emitter.emit)
        else:
            self.config: _ConfigClass = create_unbound_run_config()

        # 设置全局运行实例
        set_run(self)
        # 注册退出钩子
        self._sys_origin_excepthook = sys.excepthook
        atexit.register(self._atexit_cleanup)
        sys.excepthook = self._excepthook
        # 注册 SIGINT handler，确保 Ctrl+C 能可靠地将实验标记为 aborted
        # sys.excepthook 在主线程阻塞于 C 扩展时可能无法触发
        self._original_sigint_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self._sigint_handler)
        # 绑定日志文件，运行正式开始
        if self._ctx.config.settings.mode != "disabled":
            log.bindfile(self._ctx.debug_dir)

    # ----------------------------------
    # 私有钩子
    # ----------------------------------

    def _sigint_handler(self, signum: int, frame: Any) -> None:
        """SIGINT handler：确保 Ctrl+C 能可靠地将实验标记为 aborted。

        sys.excepthook 依赖 Python 层面抛出 KeyboardInterrupt，但当主线程
        阻塞在 C 扩展（NumPy/PyTorch 等）时，KeyboardInterrupt 可能无法
        正常传播到 excepthook。此 handler 作为额外防线，在信号层直接处理。
        """
        if self._state == "running":
            console.info("KeyboardInterrupt by user")
            import traceback

            stack = "".join(traceback.format_stack(frame)).strip() if frame is not None else ""
            error = f"KeyboardInterrupt by user\n{stack}" if stack else "KeyboardInterrupt by user"
            self.finish(state="aborted", error=error)
        # 恢复原始 handler 并重新发送信号，让进程正常终止
        signal.signal(signal.SIGINT, self._original_sigint_handler)
        if self._original_sigint_handler is signal.SIG_IGN:
            return  # Signal was ignored, do nothing.
        if callable(self._original_sigint_handler):
            self._original_sigint_handler(signum, frame)
        else:
            # The default handler (SIG_DFL) raises KeyboardInterrupt.
            raise KeyboardInterrupt

    def _atexit_cleanup(self) -> None:
        """程序正常退出时自动结束当前运行"""
        if self._state != "running":
            return
        console.debug("SwanLab Run is finishing at exit...")
        self.finish()

    def _excepthook(
        self,
        tp: Type[BaseException],
        val: BaseException,
        tb: Optional[TracebackType],
    ) -> None:
        """全局异常捕获，将实验标记为 crashed 或 aborted"""
        try:
            if self._state != "running":
                return
            state: FinishType = "crashed"
            if tp is KeyboardInterrupt:
                console.info("KeyboardInterrupt by user, aborting run...")
                state = "aborted"
            else:
                console.info("Error happened while training")
            full_error_msg = "".join(traceback.format_exception(tp, val, tb))
            self.finish(state=state, error=full_error_msg)
        except Exception as e:
            console.error(f"SwanLab failed to handle excepthook: {e}")
        finally:
            self._sys_origin_excepthook(tp, val, tb)

    def _cleanup(self):
        """
        清除副作用
        """
        # 取消钩子
        console.debug("Cleanup system hook...")
        atexit.unregister(self._atexit_cleanup)
        sys.excepthook = self._sys_origin_excepthook
        signal.signal(signal.SIGINT, self._original_sigint_handler)
        # 清理全局运行实例
        console.debug("Cleanup global instance...")
        clear_run()
        deactivate_run_config()
        console.debug("Clean & tidy! ciallo ( ∠・ω< ) ~ ★")
        # 释放日志，本次运行结束
        log.reset()

    @cached_property
    def id(self) -> str:
        """
        Current run ID.

        :return: Run ID
        """
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def run_dir(self) -> Path:
        """
        Current run directory.

        :return: Run directory path
        """
        assert self._ctx.run_dir is not None, "Run dir is not set."
        return self._ctx.run_dir

    @cached_property
    def project_url(self) -> Optional[str]:
        """
        Current project URL if in cloud mode, otherwise None.

        :return: Project URL or None
        """
        settings = self._ctx.config.settings
        if settings.mode != "cloud":
            return None
        return f"{settings.web_host}/@{settings.project.workspace}/{settings.project.name}"

    @cached_property
    def url(self) -> Optional[str]:
        """
        Current run URL if in cloud mode, otherwise None.
        :return: Run URL or None
        """
        settings = self._ctx.config.settings
        if settings.mode != "cloud":
            return None
        return f"{self.project_url}/runs/{settings.run.id}"

    # ----------------------------------
    # 上下文管理器，允许用户以 with 语句启动和结束运行
    # ----------------------------------
    def __enter__(self) -> "Run":
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if exc_type is not None:
            full_error = "".join(traceback.format_exception(exc_type, exc_val, exc_tb))
            self.finish(state="aborted" if exc_type is KeyboardInterrupt else "crashed", error=full_error)
        else:
            self.finish()

    # ----------------------------------
    # 公开 API：只负责验证输入并发事件
    # ----------------------------------
    @with_lock
    @with_run("run.log()")
    def log(self, data: Mapping[str, Any], step: Optional[int] = None):
        """Log a dictionary of metrics for the current step.

        :param data: A mapping of metric names to values. Nested dicts are flattened.

        :param step: Optional step index. If None, the step is auto-incremented.
        """
        if not (this_data := fmt.safe_validate_log_data(data)):
            console.error(f"Log data must be a dict, but got {type(data).__name__}. SwanLab will ignore this log.")
            return
        if step is not None:
            if not isinstance(step, int):
                console.error(
                    f"Step must be an integer or None, but got {type(step).__name__}. SwanLab will ignore this log."
                )
                return
            if step < 0:
                console.error(f"Step must be non-negative, but got {step}. SwanLab will ignore this log.")
                return

        next_step = self._ctx.metrics.next_step(step)

        ts = Timestamp()
        ts.GetCurrentTime()

        # 展平字典并在内部进行合规性验证和截断
        flatten_data = fmt.flatten_dict(this_data)

        # 推送日志事件
        self._emitter.emit(MetricLogEvent(data=flatten_data, step=next_step, timestamp=ts))

    @with_lock
    @with_run("run.log_scalar()")
    def log_scalar(self, *, key: str, value: Union[float, int], step: Optional[int] = None):
        """
        Log a scalar value.

        :param key: The key for the scalar value.
        :param value: The scalar value.
        :param step: Optional step for the scalar value.
        """
        self.log({key: value}, step=step)

    @with_lock
    @with_run("run.log_text()")
    def log_text(self, *, key: str, data: TextDatasType, caption: CaptionsType = None, step: Optional[int] = None):
        """
        A syntactic sugar for logging text data.

        :param key: The key for the text data.
        :param data: The text data itself or a Text object.
        :param caption: Optional caption for the text data.
        :param step: Optional step for the text data.
        """
        normalized_data = normalize_media_input(Text, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_lock
    @with_run("run.log_image()")
    def log_image(
        self,
        *,
        key: str,
        data: ImageDatasType,
        mode: ImageModesType = None,
        caption: CaptionsType = None,
        file_type: ImageFilesType = None,
        size: ImageSizesType = None,
        step: Optional[int] = None,
    ):
        """
        A syntactic sugar for logging image data.

        :param key: The key for the image data.
        :param data: The image data itself or an Image object.
        :param mode: PIL mode applied when converting to PIL.Image (e.g. 'RGB', 'L').
        :param caption: Optional caption for the image data.
        :param file_type: Output file format. One of ['png', 'jpg', 'jpeg', 'bmp']. Defaults to 'png'.
        :param size: Resize policy.
        :param step: Optional step for the image data.
        """
        normalized_data = normalize_media_input(Image, data, mode=mode, caption=caption, size=size, file_type=file_type)
        self.log({key: normalized_data}, step=step)

    @with_lock
    @with_run("run.log_audio()")
    def log_audio(
        self,
        *,
        key: str,
        data: AudioDatasType,
        sample_rate: AudioRatesType = 44100,
        caption: CaptionsType = None,
        step: Optional[int] = None,
    ):
        """
        A syntactic sugar for logging audio data.

        :param key: The key for the audio data.
        :param data: The audio data itself or an Audio object.
        :param sample_rate: Sample rate of the audio (used when data is raw numpy array).
        :param caption: Optional caption for the audio data.
        :param step: Optional step for the audio data.
        """
        normalized_data = normalize_media_input(Audio, data, caption=caption, sample_rate=sample_rate)
        self.log({key: normalized_data}, step=step)

    @with_lock
    @with_run("run.log_video()")
    def log_video(
        self,
        *,
        key: str,
        data: VideoDatasType,
        caption: CaptionsType = None,
        step: Optional[int] = None,
    ):
        """
        A syntactic sugar for logging video data.

        :param key: The key for the video data.
        :param data: The video data itself or a Video object.
        :param caption: Optional caption for the video data.
        :param step: Optional step for the video data.
        """
        normalized_data = normalize_media_input(Video, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_lock
    @with_run("run.define_scalar()")
    def define_scalar(
        self,
        *,
        key: str,
        name: Optional[str] = None,
        color: Optional[str] = None,
        x_axis: Optional[ScalarXAxisType] = None,
        chart_name: Optional[str] = None,
    ):
        """
        Manually define a scalar column before logging.

        :param key: The key for the scalar column. Supports wildcards (e.g. ``"train/*"``) to match multiple columns.
        :param name: Optional display name for the column.
        :param color: Optional color for the column, as a hex color code.
        :param x_axis: Optional x-axis type for the column.
        :param chart_name: Optional chart name to group the column into.
        """
        # TODO: 实现 glob 匹配逻辑
        if not (this_key := fmt.safe_validate_key(key)):
            return console.error(
                f"Invalid key for define scalar: {key}, please use valid characters (alphanumeric, '.', '-', '/') and avoid special characters."
            )

        original_name = name
        if name and not (name := fmt.safe_validate_name(name)):
            return console.error(f"Invalid name for define scalar: {original_name}, must be a string.")

        original_color = color
        if color and not (color := fmt.safe_validate_color(color)):
            return console.error(f"Invalid color for define scalar: {original_color}, must be a hex color code.")

        if (this_x_axis := fmt.safe_validate_x_axis(x_axis)) is None:
            return console.error(f"Invalid x_axis for define scalar: {x_axis}, must be a valid ScalarXAxisType.")

        original_chart_name = chart_name
        if chart_name and not (chart_name := fmt.safe_validate_chart_name(chart_name)):
            return console.error(f"Invalid chart_name for define scalar: {original_chart_name}, must be a string.")

        self._emitter.emit(
            ScalarDefineEvent(
                key=this_key,
                name=name,
                color=color,
                system=False,
                x_axis=this_x_axis,
                chart_name=chart_name,
                chart=None,
            )
        )

    @with_lock
    def finish(self, state: FinishType = "success", error: Optional[str] = None):
        """Finish the current run and wait for all logs to be flushed.

        :param state: Terminal state of the run. Defaults to ``"success"``.
        :param error: Optional error message, required when ``state`` is ``"crashed"``.
        """
        # 有时执行finish也有可能是系统hook主动调用，此时无需再次打印警告
        if self._state != "running":
            return
        state = state.lower()  # type: ignore
        if not (this_state := fmt.safe_validate_state(cast(FinishType, state))):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return

        self._state = this_state

        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"

        console.debug("SwanLab Run is finishing, waiting for logs to flush...")
        # 线程退出
        ts = Timestamp()
        ts.GetCurrentTime()
        self._emitter.emit(RunFinishEvent(state=this_state, error=error, timestamp=ts))
        # 阻塞主线程，等待后台队列消费完毕
        self._consumer.join()
        console.debug(f"SwanLab Run has finished with state: {self._state}, cleanup...")
        self._cleanup()


_current_run: Optional[Run] = None


def has_run() -> bool:
    """Check if there is an active SwanLab run.

    :return: True if a run is currently active, False otherwise.

    Examples:

        Check before logging:

        >>> import swanlab
        >>> if swanlab.has_run():
        ...     swanlab.log({"metric": 1.0})
        ... else:
        ...     print("No active run")
    """
    return _current_run is not None


def get_run() -> Run:
    """Get the current active SwanLab run.

    :return: The active Run instance.

    :raises RuntimeError: If no run is currently active.

    Examples:

        Access run properties:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> run = swanlab.get_run()
        >>> print(run.id)
        >>> swanlab.finish()
    """
    if _current_run is None:
        raise RuntimeError("No active Run. Call swanlab.init() first.")
    return _current_run


def set_run(run: Run) -> None:
    global _current_run
    _current_run = run


def clear_run() -> None:
    global _current_run
    _current_run = None
