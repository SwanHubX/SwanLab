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
import glob
import signal
import sys
import threading
import traceback
from concurrent.futures import Future
from functools import cached_property, wraps
from pathlib import Path
from types import TracebackType
from typing import Any, Callable, List, Literal, Mapping, Optional, Tuple, Type, Union, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    DeliverRunFinishRequest,
)
from swanlab.proto.swanlab.grpc.probe.v1.probe_pb2 import DeliverProbeStartRequest
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord
from swanlab.sdk.internal.bus import FileSaveEvent, MetricDefineEvent, MetricLogEvent
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import adapter, console, fork, helper, safe
from swanlab.sdk.internal.run import greeting
from swanlab.sdk.internal.run.components import Components
from swanlab.sdk.internal.run.components.config import Config
from swanlab.sdk.internal.run.progress import run_with_progress
from swanlab.sdk.internal.run.transforms import (
    Audio,
    ECharts,
    Html,
    Image,
    Molecule,
    Object3D,
    Text,
    Video,
    normalize_media_input,
)
from swanlab.sdk.typings.run import AsyncLogType, FinishType, ModeType, SaveType
from swanlab.sdk.typings.run.transforms import CaptionsType
from swanlab.sdk.typings.run.transforms.audio import AudioDatasType, AudioRatesType
from swanlab.sdk.typings.run.transforms.echarts import EChartsDatasType
from swanlab.sdk.typings.run.transforms.html import HtmlDatasType
from swanlab.sdk.typings.run.transforms.image import ImageDatasType, ImageFilesType, ImageModesType, ImageSizesType
from swanlab.sdk.typings.run.transforms.molecule import MoleculeDatasType
from swanlab.sdk.typings.run.transforms.object3d import Object3DDatasType
from swanlab.sdk.typings.run.transforms.text import TextDatasType
from swanlab.sdk.typings.run.transforms.video import VideoDatasType

from . import fmt

__all__ = ["Run", "has_run", "get_run", "set_run", "clear_run"]


def with_api(cmd: str, must_alive: bool = True):
    """Run API 装饰器，统一处理：fork 检测、存活校验、线程安全

    :param cmd: API 命令名，用于错误信息
    :param must_alive: True 时要求 Run 存活，False 时仅做线程安全（如 finish 允许在非存活状态调用）
    """

    def decorator(f):
        @wraps(f)
        def wrapper(self: "Run", *args, **kwargs):
            if fork.is_forked(self._init_pid):
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

        url: Online URL for this run (None in local/offline/disabled mode).

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

    def __init__(self, ctx: RunContext, path: Optional[str] = None):
        # 1. 基础状态、组件准备
        self._ctx = ctx
        self._is_main_thread = threading.current_thread() is threading.main_thread()
        # 运行实验路径，/:username/:project_name/:run_id
        self._path = path
        self._state: Union[FinishType, Literal["running"]] = "running"
        self._init_pid = fork.current_pid()
        # 外部API锁，防止并发调用
        self._api_lock = threading.RLock()
        # 运行时组件
        self._components = Components(ctx)
        # 回调器
        self._callbacker = self._ctx.callbacker
        # 独立组件
        self._core = self._ctx.core
        self._probe = self._ctx.probe
        # 2. 注册副作用
        # 设置全局运行实例
        set_run(self)
        # TODO: swanlab-core 上线后，此处注册 Run 级别的 fork 重连回调
        # 注册退出钩子
        self._sys_origin_excepthook = sys.excepthook
        atexit.register(self._handle_atexit)
        sys.excepthook = self._handle_except
        # 注册 SIGINT handler，确保 Ctrl+C 能可靠地将实验标记为 aborted，sys.excepthook 在主线程阻塞于 C 扩展时可能无法触发
        # signal() 仅允许在主线程调用；Ray actor 等非主线程环境下跳过注册，
        # atexit / excepthook 钩子仍然生效
        self._original_sigint_handler = signal.getsignal(signal.SIGINT)
        if self._is_main_thread:
            signal.signal(signal.SIGINT, self._handle_sigint)

        # 3. 启动组件 + 初始化日志
        # 回调的path在path为空的时候自动生成一个/:project_name/:run_id，否则使用path
        run_settings = ctx.config.settings
        assert run_settings.project.name, "project name is required when init run"
        assert run_settings.run.id, "run id is required when init run"
        run_path = path or f"/{run_settings.project.name}/{run_settings.run.id}"
        callback_settings = run_settings.model_copy(deep=True)
        self._callbacker.on_run_initialized(self._ctx.run_dir, run_path, settings=callback_settings)
        # 启动组件
        self._components.start()
        # 启动硬件监控探针
        start_request = DeliverProbeStartRequest(
            probe_settings=run_settings.to_probe_proto(
                run_id=run_settings.run.id,
                run_dir=self._ctx.run_dir,
                global_system_step=self._ctx.global_system_step,
            )
        )
        self._probe.deliver_probe_start(start_request)
        console.init(bind_to=self._ctx.debug_dir if self.mode != "disabled" else None)
        greeting.welcome(self._ctx, self)

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
        Current run path in the format of /:username/:project_name/:run_id.

        :return: Run path
        """
        if not self._path:
            raise ValueError(
                "Run path is unavailable because the current run is not using SwanLab online mode. "
                "Please initialize SwanLab with online mode to access the run path."
            )
        return self._path

    @cached_property
    def url(self) -> str:
        """
        Current run URL if in online mode, otherwise None.
        :return: Run URL or None
        """
        if not self._path:
            raise ValueError(
                "Run url is unavailable because the current run is not using SwanLab online mode. "
                "Please initialize SwanLab with online mode to access the run path."
            )
        settings = self._ctx.config.settings
        return f"{settings.web_host}{helper.fmt_run_path(self._path)}"

    @cached_property
    def config(self) -> Config:
        return self._components.config

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
        return not fork.is_forked(self._init_pid) and self._state == "running"

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

    @safe.decorator(message="SwanLab failed to log data due to unexpected error")
    def _log_impl(self, data: Mapping[str, Any], step: Optional[int] = None):
        """log 的无锁内部实现，供 async_log 回调调用以避免与 finish() 的 _api_lock 死锁。"""
        # 1. 验证输入
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

        # 2. 获取当前log快照，如步骤、时间戳等
        next_step = self._ctx.next_step(step)
        ts = Timestamp()
        ts.GetCurrentTime()

        # 3. 推送日志事件
        # 展平字典并在内部进行合规性验证和截断
        flatten_data = fmt.flatten_dict(this_data)
        self._components.emitter.emit(MetricLogEvent(data=flatten_data, step=next_step, timestamp=ts))

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

        return self._components.asynctask.submit(
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

    @with_api("run.log_echarts()")
    def log_echarts(
        self, *, key: str, data: EChartsDatasType, caption: CaptionsType = None, step: Optional[int] = None
    ):
        """
        A syntactic sugar for logging ECharts data.

        :param key: The key for the echarts data.
        :param data: The echarts chart object (pyecharts chart) or an ECharts object.
        :param caption: Optional caption for the echarts data.
        :param step: Optional step for the echarts data.
        """
        normalized_data = normalize_media_input(ECharts, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_api("run.log_object3d()")
    def log_object3d(
        self, *, key: str, data: Object3DDatasType, caption: CaptionsType = None, step: Optional[int] = None
    ):
        """
        A syntactic sugar for logging 3D object data.

        :param key: The key for the 3D object data.
        :param data: The 3D object data (numpy array, dict, file path, or Object3D object).
        :param caption: Optional caption for the 3D object data.
        :param step: Optional step for the 3D object data.
        """
        normalized_data = normalize_media_input(Object3D, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_api("run.log_molecule()")
    def log_molecule(
        self, *, key: str, data: MoleculeDatasType, caption: CaptionsType = None, step: Optional[int] = None
    ):
        """
        A syntactic sugar for logging molecule data.

        :param key: The key for the molecule data.
        :param data: The molecule data (SMILES string, file path, RDKit Mol object, or Molecule object).
        :param caption: Optional caption for the molecule data.
        :param step: Optional step for the molecule data.
        """
        normalized_data = normalize_media_input(Molecule, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_api("run.log_html()")
    def log_html(self, *, key: str, data: HtmlDatasType, caption: CaptionsType = None, step: Optional[int] = None):
        """
        A syntactic sugar for logging HTML data.

        :param key: The key for the HTML data.
        :param data: The HTML data (raw HTML string, file path, file object, or Html object).
        :param caption: Optional caption for the HTML data.
        :param step: Optional step for the HTML data.
        """
        normalized_data = normalize_media_input(Html, data, caption=caption)
        self.log({key: normalized_data}, step=step)

    @with_api("run.define_metric()")
    def define_metric(
        self,
        key: str,
        *,
        x_axis: Optional[str] = None,
        section_name: Optional[str] = None,
        hidden: Optional[bool] = None,
        step_sync: Optional[bool] = None,
        overwrite: bool = False,
        **kwargs,
    ):
        """Define a metric's display configuration before logging.

        Customizes how an auto-generated chart for ``key`` appears in project Views:
        X-axis, section placement, and visibility. The same ``(class, key)``
        shares one chart across all runs in the project.

        :param key: Metric key. Supports exact match and a single trailing ``*``
            glob (e.g. ``"train/*"``). System keys are never matched.
        :param x_axis: Custom X-axis key. ``None`` (default) means the system
            step. The system axes ``"_step"`` and ``"_relative_time"`` are also
            accepted: they are resolved by the server, and the SDK performs no
            step injection or X dedup for them. Only affects scalar charts;
            media ignores this. The X series is assumed monotonically
            non-decreasing — for a given X value, only the first logged Y is
            kept (consecutive-duplicate X values are dropped). ``step_metric``
            is accepted as a ``**kwargs`` alias for backward compatibility.
        :param section_name: Section name for the auto chart. ``None`` means
            the default section derived from the key.
        :param hidden: If ``True``, place the chart in the HIDDEN section.
            Three states: ``None`` (default) means "not provided" — merge mode
            keeps the previous value; ``True`` hides; ``False`` explicitly
            unhides (also effective in merge mode).
        :param step_sync: Whether this Y key should copy the latest custom X
            value onto the current step when X and Y are logged separately.
            Defaults to ``True`` when ``x_axis`` (or ``step_metric``) is set,
            ``False`` otherwise — without a custom X axis it has no effect.
            ``False`` means this Y will not trigger X injection: the two series
            only align when they share a step (same ``log()`` or the same
            explicit ``step``). Duplicate-X dropping then uses only an X value
            present in this event. Sibling Y keys that keep ``step_sync=True``
            can still inject X for the shared series.
        :param overwrite: If ``False`` (default), merge with previous calls for
            the same ``key`` — unspecified fields reuse the previous value.
            If ``True``, unspecified fields reset to their default, overwriting
            previous values. Only affects rules not yet applied to a logged key.

        .. note::
            Project-wide first-definition-effective. A chart is shared across every run
            in the project under the same ``(class, key)``, and
            ``define_metric`` only shapes it on the **first** run that
            introduces the key to the project. Once the chart exists — created
            by this run's first log or by any earlier run — later
            ``define_metric`` calls (including ``overwrite=True``) have no
            effect, and no warning is raised. To change an existing chart,
            edit it in the UI.

        Examples:

            Custom X-axis with separate logging:

            >>> import swanlab
            >>> swanlab.init(mode="local")
            >>> run = swanlab.get_run()
            >>> run.define_metric("train/loss", x_axis="train/epoch")
            >>> swanlab.log({"train/epoch": 1})    # step=1
            >>> swanlab.log({"train/loss": 0.8})   # step=2, auto-syncs train/epoch=1

            Glob pattern for all validation metrics:

            >>> run.define_metric("val/*", section_name="Validation")

            Disable X injection (X/Y must share a step to align):

            >>> run.define_metric("train/acc", x_axis="train/epoch", step_sync=False)
            >>> swanlab.log({"train/epoch": 1}, step=1)
            >>> swanlab.log({"train/acc": 0.9}, step=2)  # no injected epoch at step=2

            ``step_metric`` alias:

            >>> run.define_metric("loss", step_metric="custom_step")
        """
        # 1. key 校验，完全复用 log 侧的 validate_key（非字符串强转 str、清洗首尾空白与'./'、超长截断）
        # 保证规则与 log 实际产生的规范 key 落在同一形式，避免发出永不匹配的规则
        #
        # 边界情况：>512 的 key 截断可能截掉末尾 glob '*'，glob 退化为 exact 匹配，此时属病态输入且 log 侧同样截断、两边落点一致，暂不特殊处理
        try:
            this_key = fmt.validate_key(key)
        except ValueError as e:
            console.error(f"Invalid key for define_metric: {key!r}, {e}")
            return

        # 2. glob 校验与分类：仅支持末尾单个 '*'（如 train/*），拒绝 *loss、train/*/loss、train/** 等
        star_count = this_key.count("*")
        if star_count > 0 and not (star_count == 1 and this_key.endswith("*")):
            console.error(
                f"Invalid glob pattern for define_metric: {key!r}, "
                f"only a single trailing '*' is supported (e.g. 'train/*')"
            )
            return
        is_glob = star_count == 1

        # 3. step_metric 兼容别名（仅从 kwargs 读取；其余未知 kwargs 静默忽略，与 init 兼容风格一致）
        step_metric = kwargs.pop("step_metric", None)
        if step_metric is not None:
            if x_axis is None:
                x_axis = step_metric
            elif x_axis != step_metric:
                console.warning(
                    f"Conflicting x_axis={x_axis!r} and step_metric={step_metric!r} "
                    f"for key {this_key!r}, using x_axis={x_axis!r}"
                )

        # 4. x_axis / section_name 校验
        if x_axis is not None:
            validated_x = fmt.safe_validate_x_axis(x_axis)
            if validated_x is None:
                console.error(
                    f"Invalid x_axis for define_metric: {x_axis}, must be a valid metric key, "
                    f"'_step', '_relative_time', and must not be a system metric key."
                )
                return
            x_axis = validated_x
        # section_name：非法 / 空串降级为默认 section（不阻断整条 define）
        if section_name is not None:
            validated_section = fmt.safe_validate_name(section_name)
            if validated_section is None or not validated_section.strip():
                console.warning(
                    f"Invalid section_name for define_metric: {section_name!r}, ignored; using default section."
                )
                section_name = None
            else:
                section_name = validated_section

        # 5. 发射事件（三态字段校验后即为 None|str / None|bool，直接透传；is_glob 已在第 2 步分类）
        self._components.emitter.emit(
            MetricDefineEvent(
                key=this_key,
                is_glob=is_glob,
                x_axis=x_axis,
                section_name=section_name,
                hidden=hidden,
                step_sync=step_sync,
                overwrite=overwrite,
            )
        )

    @with_api("run.save()")
    def save(
        self,
        glob_str: Union[str, bytes, Path],
        base_path: Optional[Union[str, Path]] = None,
        policy: Literal["now", "end", "live"] = "live",
    ) -> List[str]:
        """Save files matched by glob into the current run.

        :param glob_str: A glob pattern matching files to save (e.g. ``"checkpoints/*.pt"``).
        :param base_path: Base directory for resolving relative paths. Defaults to cwd.
        :param policy: Save policy:

            - ``"now"`` — upload matched files immediately.
            - ``"end"`` — defer upload until the run finishes.
            - ``"live"`` — watch for file changes and re-upload automatically.

        :return: List of matched file paths (relative to base_path).
        """
        if not (this_policy := fmt.safe_validate_save_policy(policy)):
            console.error(f"Invalid save policy: {policy}, allowed values are {get_args(SaveType)}")
            return []

        resolved_paths = fmt.resolve_save_paths(glob_str, base_path)
        if resolved_paths is None:
            return []
        resolved_glob, resolved_base = resolved_paths

        # Glob 匹配文件
        matched = glob.glob(str(resolved_glob), recursive=True)
        if not matched:
            console.warning(f"No files matched by glob pattern: {glob_str}")
            return []

        # 过滤出普通文件，计算相对路径
        core_settings = self._ctx.config.settings.core
        files: List[Tuple[Path, Path]] = []
        for abs_str in matched:
            abs_path = Path(abs_str)
            if not abs_path.is_file():
                continue
            try:
                rel_path = abs_path.relative_to(resolved_base)
            except ValueError:
                continue
            files.append((abs_path, rel_path))

        if not files:
            console.warning(f"No files matched by glob pattern: {glob_str}")
            return []

        # 校验批次大小
        if len(files) > core_settings.save_batch:
            raise ValueError(f"Too many files matched ({len(files)}), limit is {core_settings.save_batch}")

        # 按 policy 分发事件
        # Core 负责创建 symlink 处理和 end policy 下的延迟上传
        results: List[str] = []
        for source_path, rel_path in files:
            event = FileSaveEvent(
                source_path=str(source_path),
                name=str(rel_path),
                policy=this_policy,
            )
            self._components.emitter.emit(event)
            results.append(str(rel_path))
        return results

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
        if not self.alive:
            console.warning("SwanLab Run has already finished or has not started.")
            return
        state = state.lower()  # type: ignore
        if not (this_state := fmt.safe_validate_state(cast(FinishType, state))):
            console.error(f"Invalid state: {state}, allowed values are {get_args(FinishType)}")
            return
        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"
        greeting.goodbye(self._ctx, self)
        # 2. 运行结束，清理组件
        self._state = this_state
        # 停止所有内部组件（async_log → terminal → config → consumer）
        self._components.stop(async_log_timeout=async_log_timeout)
        # 停止probe
        self._probe.deliver_probe_finish()
        ts = Timestamp()
        ts.GetCurrentTime()
        # 3. 停止Core
        finish_request = DeliverRunFinishRequest(
            finish_record=FinishRecord(state=adapter.state[this_state], error=error, finished_at=ts)
        )
        finish_resp = self._core.deliver_run_finish(finish_request)
        if not finish_resp.success:
            console.error(finish_resp.message)
        confirm_resp = run_with_progress(
            stats_fn=self._core.get_operation_stats,
            blocking_fn=self._core.confirm_run_finish,
        )
        if not confirm_resp.success:
            console.error(confirm_resp.message)
        # finish 回调
        self._callbacker.on_run_finished(this_state, error)
        console.debug(f"SwanLab Run has finished with state: {self._state}, cleanup...")
        # 4. 清理副作用
        console.debug("Cleanup system hook...")
        atexit.unregister(self._handle_atexit)
        sys.excepthook = self._sys_origin_excepthook
        if self._is_main_thread:
            signal.signal(signal.SIGINT, self._original_sigint_handler)
        # 清理全局运行实例
        console.debug("Cleanup global instance...")
        clear_run()
        # [随临时方案删除] online 模式下销毁全局 client 单例，确保下次 init() 重新认证获取新 sid，
        # 避免服务端在实验结束后使 sid 失效导致 401；待 client 生命周期归属 Core 后删除，跟踪 issue: #1742
        if self.mode == "online" and client.exists():
            console.debug("Reset online client...")
            client.reset()
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
