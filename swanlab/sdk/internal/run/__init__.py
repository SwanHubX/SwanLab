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
import os
import signal
import sys
import threading
import traceback
from concurrent.futures import Future
from functools import cached_property, wraps
from pathlib import Path
from types import TracebackType
from typing import Any, Callable, Literal, Mapping, Optional, Type, Union, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord
from swanlab.sdk.internal.bus.events import (
    MetricLogEvent,
    ScalarDefineEvent,
)
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import adapter, console, safe
from swanlab.sdk.internal.run import system
from swanlab.sdk.internal.run.async_task import AsyncTaskManager
from swanlab.sdk.internal.run.config import (
    Config,
    deactivate_run_config,
)
from swanlab.sdk.internal.run.transforms import Audio, Image, Text, Video, normalize_media_input
from swanlab.sdk.typings.run import AsyncLogType, FinishType, ModeType
from swanlab.sdk.typings.run.column import ScalarXAxisType
from swanlab.sdk.typings.run.transforms import CaptionsType
from swanlab.sdk.typings.run.transforms.audio import AudioDatasType, AudioRatesType
from swanlab.sdk.typings.run.transforms.image import ImageDatasType, ImageFilesType, ImageModesType, ImageSizesType
from swanlab.sdk.typings.run.transforms.text import TextDatasType
from swanlab.sdk.typings.run.transforms.video import VideoDatasType

from . import utils_fmt as fmt
from .factory import factory_config, factory_consumer, factory_emitter, factory_monitor
from .record_builder import RecordBuilder

__all__ = ["Run", "has_run", "get_run", "set_run", "clear_run"]


def with_api(cmd: str, must_alive: bool = True):
    """Run API 装饰器，统一处理：fork 检测、存活校验、线程安全

    :param cmd: API 命令名，用于错误信息
    :param must_alive: True 时要求 Run 存活，False 时仅做线程安全（如 finish 允许在非存活状态调用）
    """

    def decorator(f):
        @wraps(f)
        def wrapper(self: "Run", *args, **kwargs):
            if self._forked:
                # fork 后子进程继承的锁可能已持有，替换为新锁避免死锁
                self._api_lock = threading.RLock()
                raise RuntimeError(
                    "SwanLab Run does not support fork yet. Use `multiprocessing.set_start_method('spawn')` "
                    "or call `swanlab.init()` in the child process."
                )
            if must_alive and not self.alive:
                raise RuntimeError(f"`{cmd}` requires an active Run, call `swanlab.init()` first.")
            with self._api_lock:
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

        dir: Local directory where run data is stored.

        url: Cloud URL for this run (None in local/offline/disabled mode).

    Examples:

        Access run properties:

        >>> import swanlab
        >>> run = swanlab.init(mode="local", project="my_project")
        >>> print(run.id)
        >>> print(run.dir)
        >>> swanlab.finish()

        Log metrics:

        >>> import swanlab
        >>> run = swanlab.init(mode="local")
        >>> run.log({"loss": 0.5, "accuracy": 0.95})
        >>> swanlab.finish()
    """

    def __init__(self, ctx: RunContext):
        # 1. 基础状态、组件准备
        self._ctx = ctx
        self._state: Union[FinishType, Literal["running"]] = "running"
        self._pid = os.getpid()
        # 外部API锁，防止并发调用
        self._api_lock = threading.RLock()
        # 异步任务管理器：处理async_log任务
        self._async_task_manager = AsyncTaskManager()
        # 运行时组件
        self._builder = RecordBuilder(self._ctx)
        self._emitter = factory_emitter(self._ctx)
        self._config = factory_config(self._ctx, self._emitter)
        self._consumer = factory_consumer(self._ctx, self._emitter, self._builder)
        self._callbacker = self._ctx.callbacker
        # self._monitor is not None 则代表硬件监控开启
        self._monitor: Optional[system.Monitor] = None

        # 2. 注册副作用
        # 设置全局运行实例
        set_run(self)
        # 注册退出钩子
        self._sys_origin_excepthook = sys.excepthook
        atexit.register(self._handle_atexit)
        sys.excepthook = self._handle_except
        # 注册 SIGINT handler，确保 Ctrl+C 能可靠地将实验标记为 aborted，sys.excepthook 在主线程阻塞于 C 扩展时可能无法触发
        self._original_sigint_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self._handle_sigint)

        # 3. 初始化完成
        self._callbacker.on_run_initialized(self._ctx.run_dir, self.path)
        self._monitor = factory_monitor(self._ctx, self._emitter)
        self._consumer.start()
        # 初始化日志模块：非 disabled 模式绑定文件，disabled 模式禁用持久化
        console.init(bind_to=self._ctx.debug_dir if self.mode != "disabled" else None)

    # ----------------------------------
    # 私有钩子
    # ----------------------------------

    def _handle_sigint(self, signum: int, frame: Any) -> None:
        """SIGINT handler：确保 Ctrl+C 能可靠地将实验标记为 aborted。

        sys.excepthook 依赖 Python 层面抛出 KeyboardInterrupt，但当主线程阻塞在 C 扩展（NumPy/PyTorch 等）时，
        KeyboardInterrupt 可能无法正常传播到 excepthook。
        此 handler 作为额外防线，在信号层直接处理。
        """
        if self.alive:
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

    def _handle_atexit(self) -> None:
        """程序正常退出时自动结束当前运行"""
        if not self.alive:
            return
        console.debug("SwanLab Run is finishing at exit...")
        self.finish()

    def _handle_except(self, tp: Type[BaseException], val: BaseException, tb: Optional[TracebackType]) -> None:
        """全局异常捕获，将实验标记为 crashed 或 aborted"""
        with safe.block(message="SwanLab failed to handle excepthook"):
            if self.alive:
                state: FinishType = "crashed"
                if tp is KeyboardInterrupt:
                    console.info("KeyboardInterrupt by user, aborting run...")
                    state = "aborted"
                else:
                    console.info("Error happened while training")
                full_error_msg = "".join(traceback.format_exception(tp, val, tb))
                self.finish(state=state, error=full_error_msg)
        self._sys_origin_excepthook(tp, val, tb)

    # ----------------------------------
    # 公开辅助属性
    # ----------------------------------

    @cached_property
    def id(self) -> str:
        """
        Current run ID.

        :return: Run ID
        """
        assert self._ctx.config.settings.run.id is not None, "Run id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def mode(self) -> ModeType:
        """
        Current run mode.

        :return: Run mode
        """
        assert self._ctx.config.settings.mode is not None, "Run mode is not set."
        return self._ctx.config.settings.mode

    @cached_property
    def name(self) -> str:
        """
        Current run name, equal to experiment name.
        """
        assert self._ctx.config.settings.experiment.name is not None, "Experiment name is not set."
        return self._ctx.config.settings.experiment.name

    @cached_property
    def dir(self) -> Path:
        """
        Current run directory.

        :return: Run directory path
        """
        assert self._ctx.run_dir is not None, "Run dir is not set."
        return self._ctx.run_dir

    @cached_property
    def path(self) -> str:
        """
        Current run path in the format of /:project/:run_id.

        :return: Run path
        """
        settings = self._ctx.config.settings
        return f"/{settings.project.name}/runs/{settings.run.id}"

    @cached_property
    def url(self) -> Optional[str]:
        """
        Current run URL if in cloud mode, otherwise None.
        :return: Run URL or None
        """
        settings = self._ctx.config.settings
        if settings.mode != "cloud":
            return None
        return f"{settings.web_host}/@{settings.project.workspace}{self.path}"

    @cached_property
    def config(self) -> Config:
        return self._config

    @property
    def _forked(self) -> bool:
        """当前进程是否为创建 Run 时的进程的 fork 子进程"""
        return os.getpid() != self._pid

    @property
    def _passive(self) -> bool:
        """被动模式：仅解析验证，不产生任何副作用（文件IO、网络、线程等）。

        当前 disabled 模式为被动模式，未来其他需要跳过运行时组件的模式也可复用此属性。
        """
        return self.mode == "disabled"

    @property
    def alive(self) -> bool:
        """
        If the run is alive. You can log metrics if the run is alive.
        :return: True if the run is alive, False otherwise
        """
        return not self._forked and self._state == "running"

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
    @with_api("run.log()")
    def log(self, data: Mapping[str, Any], step: Optional[int] = None):
        """Log a dictionary of metrics for the current step.

        :param data: A mapping of metric names to values. Nested dicts are flattened.

        :param step: Optional step index. If None, the step is auto-incremented.
        """
        self._log_impl(data, step)

    def _log_impl(self, data: Mapping[str, Any], step: Optional[int] = None):
        """log 的无锁内部实现，供 async_log 回调调用以避免与 finish() 的 _api_lock 死锁。"""
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

    @with_api("run.async_log()")
    def async_log(
        self,
        func: Callable,
        *args,
        step: Optional[int] = None,
        mode: AsyncLogType = "threading",
        **kwargs,
    ) -> Future:
        """Asynchronously execute a function and automatically log its return value.

        ``func`` is submitted to a background thread, process, or the asyncio event loop (depending on *mode*).
        When it completes, its return value — a ``dict`` — is passed to :meth:`log` automatically.
        The call returns a :class:`~concurrent.futures.Future` immediately.

        ``finish()`` waits for all outstanding ``async_log`` tasks before flushing, so no data is lost.

        :param func: A callable returning a ``dict`` suitable for :meth:`log`.
        :param args: Positional arguments forwarded to *func*.
        :param step: Optional step index. If ``None``, auto-incremented when the task **completes** (not when
            submitted). Pass an explicit value if step ordering matters.
        :param mode: Execution mode:

            - ``"asyncio"`` — schedule on the running asyncio event loop. *func* must be a coroutine
              (``async def``). No pickle constraints. Raises :exc:`RuntimeError` if no loop is running.

            - ``"threading"`` (default) — background thread. No pickle constraints; *func* can access
              ``swanlab.config`` and return media objects (:class:`Image`, :class:`Audio`, etc.).
              Subject to the GIL.

            - ``"spawn"`` — new child process (``mp_context=spawn``). Bypasses the GIL, ideal for CPU-bound
              work. *func*, its arguments, and its return value **must be pickle-serializable** (no
              :class:`Image`, ``torch.Tensor``, etc.). The child process cannot access the active Run.

            - ``"fork"`` — **reserved**. Will be enabled after ``swanlab-core`` ships; forked children
              will call ``swanlab.log()`` directly, removing the pickle constraint.

        :param kwargs: Keyword arguments forwarded to *func*.
        :return: A :class:`~concurrent.futures.Future`. In ``"asyncio"`` mode the future is asyncio-compatible
            (wrapped via :func:`asyncio.wrap_future`).
        :raises RuntimeError: No active Run, or no asyncio event loop (``"asyncio"`` mode only).

        Examples:

            Asyncio mode — coroutine function for IO-bound work:

            >>> import swanlab
            >>> run = swanlab.init()
            >>> async def slow_compute():
            ...     import asyncio
            ...     await asyncio.sleep(2)
            ...     return {"score": 0.95}
            >>> future = run.async_log(slow_compute, step=1, mode="asyncio")

            Threading mode (default) — IO-bound or returning media objects:

            >>> def fetch_score():
            ...     import time, numpy as np
            ...     time.sleep(2)
            ...     return {"score": 0.95, "preview": swanlab.Image(np.random.randn(10, 10))}
            >>> future = run.async_log(fetch_score, step=1)

            Spawn mode — CPU-bound, pickle-safe return values:

            >>> def compute_loss():
            ...     return {"loss": 0.123, "acc": 0.95}
            >>> future = run.async_log(compute_loss, step=2, mode="spawn")

            Spawn mode with torch — convert before returning:

            >>> def compute():
            ...     import torch
            ...     t = torch.randn(10)
            ...     return {"value": t.item(), "arr": t.detach().cpu().numpy()}
            >>> future = run.async_log(compute, step=3, mode="spawn")
        """
        if mode == "fork":
            raise NotImplementedError(
                "fork mode is not yet supported, please looking forward to the `swanlab-core` release"
            )

        return self._async_task_manager.submit(
            func,
            args=args,
            kwargs=kwargs,
            step=step,
            mode=mode,
            on_success=lambda result, s: self._log_impl(result, step=s),
            on_error=lambda: console.trace("swanlab.async_log run error"),
        )

    @with_api("run.log_scalar()")
    def log_scalar(self, *, key: str, value: Union[float, int], step: Optional[int] = None):
        """
        Log a scalar value.

        :param key: The key for the scalar value.
        :param value: The scalar value.
        :param step: Optional step for the scalar value.
        """
        self.log({key: value}, step=step)

    @with_api("run.log_text()")
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

    @with_api("run.log_image()")
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

    @with_api("run.log_audio()")
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

    @with_api("run.log_video()")
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

    @with_api("run.define_scalar()")
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

    @with_api("run.finish()", must_alive=False)
    def finish(
        self,
        state: FinishType = "success",
        error: Optional[str] = None,
        async_log_timeout: Optional[int] = None,
    ):
        """Finish the current run and wait for all logs to be flushed.

        :param state: Terminal state of the run. Defaults to ``"success"``.
        :param error: Optional error message, required when ``state`` is ``"crashed"``.
        :param async_log_timeout: Optional timeout for async_log tasks. None means no timeout.
        """
        # 1. 状态校验
        # 有时执行finish也有可能是系统hook主动调用，此时无需再次打印警告，如果在finishing状态，也忽略
        if not self.alive:
            return
        state = state.lower()  # type: ignore
        if not (this_state := fmt.safe_validate_state(cast(FinishType, state))):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return
        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"
        # 2. 运行结束前，结束其他依赖于运行实例的线程
        # 2.1 等待所有 async_log 任务完成
        console.debug("Waiting for async_log tasks to complete...")
        self._async_task_manager.shutdown(timeout=async_log_timeout)
        # 2.2 停止硬件监控
        if self._monitor is not None:
            console.debug("Stopping hardware monitor...")
            self._monitor.stop()

        # 3. 运行结束
        self._state = this_state
        # 停止时间
        ts = Timestamp()
        ts.GetCurrentTime()
        # 3.1 TODO: goodbye message

        # 3.2 TODO: 停止终端代理

        # 3.3 停止消费者线程
        console.debug("SwanLab Run is finishing, waiting for logs to flush...")
        self._consumer.stop()
        self._consumer.join()
        # 3.4 停止Core线程
        finish_resp = self._ctx.core.deliver_run_finish(
            FinishRecord(state=adapter.state[this_state], error=error, finished_at=ts)
        )
        if not finish_resp.success:
            console.error(finish_resp.message)
        console.debug(f"SwanLab Run has finished with state: {self._state}, cleanup...")
        # 3.5 清理副作用
        console.debug("Cleanup system hook...")
        atexit.unregister(self._handle_atexit)
        sys.excepthook = self._sys_origin_excepthook
        signal.signal(signal.SIGINT, self._original_sigint_handler)
        # 清理全局运行实例
        console.debug("Cleanup global instance...")
        clear_run()
        deactivate_run_config()
        console.debug("Clean & tidy! ciallo ( ∠・ω< ) ~ ★")
        # 释放日志，本次运行结束
        console.reset()


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
    return _current_run is not None and _current_run.alive


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
