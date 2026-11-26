"""
@author: mx
@file: eval.py
@time: 2026/11/25
@description: 评估运行时对象与全局槽位，管理一次 swanlab.eval() 的完整生命周期
"""

from __future__ import annotations

import atexit
import signal
import sys
import threading
import traceback
from collections.abc import Mapping
from functools import cached_property
from pathlib import Path
from types import TracebackType
from typing import Any, Dict, Literal, Optional, Type, cast, get_args

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.eval.v1.eval_pb2 import EvalFinishRecord
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import DeliverEvalFinishRequest
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import adapter, console, fork, safe
from swanlab.sdk.typings.run import FinishType, ModeType

from .benchmark import MAX_BENCHMARK_NAME_LENGTH, Benchmark
from .components import Components
from .components.config import Config, parse
from .progress import run_with_progress

__all__ = ["EvalRun", "has_eval_run", "get_eval_run", "set_eval_run", "clear_eval_run"]

# 合法的结束状态，与 FinishType 保持一致
_FINISH_STATES = set(get_args(FinishType))


class EvalRun:
    """A SwanLab evaluation run.

    由 ``swanlab.eval()`` 创建，持有本次评估的 Core、配置与评测项注册表。
    通过 ``benchmark()`` 声明评测项并调用 ``log()`` 上报分数。

    Attributes:

        config: Configuration object for storing evaluation hyperparameters.

        id: Unique identifier of this eval.

        name: Eval name.

        dir: Local directory where eval data is stored.

    Examples:

        >>> import swanlab
        >>> with swanlab.eval(project="model-eval", name="checkpoint-5000") as run:
        ...     mmlu = run.benchmark("MMLU", config={"num_fewshot": 5})
        ...     mmlu.log({"accuracy": 0.72})
    """

    def __init__(self, ctx: RunContext):
        self._ctx = ctx
        self._state: Literal["running", "finished"] = "running"
        self._init_pid = fork.current_pid()
        self._is_main_thread = threading.current_thread() is threading.main_thread()
        # 评测项注册与运行结束互斥
        self._api_lock = threading.RLock()
        self._benchmarks: Dict[str, Benchmark] = {}
        self._components = Components(ctx)
        self._core = ctx.core
        # 设置全局评估实例
        set_eval_run(self)
        # 注册退出钩子
        self._sys_origin_excepthook = sys.excepthook
        atexit.register(self._handle_atexit)
        sys.excepthook = self._handle_except
        # signal() 仅允许在主线程调用
        self._original_sigint_handler = signal.getsignal(signal.SIGINT)
        if self._is_main_thread:
            signal.signal(signal.SIGINT, self._handle_sigint)
        # 启动组件 + 初始化日志
        self._components.start()
        run_settings = ctx.config.settings
        assert run_settings.project.name, "project name is required when eval starts"
        assert run_settings.run.id, "eval id is required when eval starts"
        if self.mode == "disabled" or run_settings.core.skip_store:
            console.init(bind_to=None)
        else:
            console.init(bind_to=self._ctx.debug_dir)

    # ----------------------------------
    # 私有钩子
    # ----------------------------------

    def _handle_atexit(self) -> None:
        """程序正常退出时自动结束当前评估"""
        if not self.alive:
            return
        console.debug("SwanLab Eval is finishing at exit...")
        self.finish()

    def _handle_sigint(self, signum: int, frame: Any) -> None:
        """SIGINT handler：确保 Ctrl+C 能记录为中断。"""
        if self.alive:
            console.info("KeyboardInterrupt by user")
            stack = "".join(traceback.format_stack(frame)).strip() if frame is not None else ""
            error = f"KeyboardInterrupt by user\n{stack}" if stack else "KeyboardInterrupt by user"
            self.finish(state="aborted", error=error)
        # 恢复原始 handler 并重新发送信号，让进程正常终止
        signal.signal(signal.SIGINT, self._original_sigint_handler)
        if self._original_sigint_handler is signal.SIG_IGN:
            return
        if callable(self._original_sigint_handler):
            self._original_sigint_handler(signum, frame)
        else:
            raise KeyboardInterrupt

    def _handle_except(self, tp: Type[BaseException], val: BaseException, tb: Optional[TracebackType]) -> None:
        """全局异常捕获，将评估标记为 crashed 或 aborted"""
        with safe.block(message="SwanLab failed to handle eval excepthook"):
            if self.alive:
                state: FinishType = "crashed"
                if tp is KeyboardInterrupt:
                    console.info("KeyboardInterrupt by user, aborting eval...")
                    state = "aborted"
                else:
                    console.info("Error happened while evaluating")
                full_error_msg = "".join(traceback.format_exception(tp, val, tb))
                self.finish(state=state, error=full_error_msg)
        self._sys_origin_excepthook(tp, val, tb)

    # ----------------------------------
    # 公开辅助属性
    # ----------------------------------

    @cached_property
    def id(self) -> str:
        """
        Current eval id.

        :return: Eval id
        """
        assert self._ctx.config.settings.run.id is not None, "Eval id is not set."
        return self._ctx.config.settings.run.id

    @cached_property
    def mode(self) -> ModeType:
        """
        Current eval mode.

        :return: Eval mode
        """
        assert self._ctx.config.settings.mode is not None, "Eval mode is not set."
        return self._ctx.config.settings.mode

    @cached_property
    def name(self) -> str:
        """
        Current eval name.
        """
        assert self._ctx.config.settings.experiment.name is not None, "Eval name is not set."
        return self._ctx.config.settings.experiment.name

    @cached_property
    def dir(self) -> Path:
        """
        Current eval directory.

        :return: Eval directory path
        """
        assert self._ctx.run_dir is not None, "Eval dir is not set."
        return self._ctx.run_dir

    @cached_property
    def config(self) -> Config:
        return self._components.config

    @property
    def alive(self) -> bool:
        """
        If the eval is alive. You can log scores if the eval is alive.
        :return: True if the eval is alive, False otherwise
        """
        return not fork.is_forked(self._init_pid) and self._state == "running"

    # ----------------------------------
    # 上下文管理器，允许用户以 with 语句启动和结束评估
    # ----------------------------------
    def __enter__(self) -> "EvalRun":
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if exc_type is not None:
            full_error = "".join(traceback.format_exception(exc_type, exc_val, exc_tb))
            self.finish(state="aborted" if exc_type is KeyboardInterrupt else "crashed", error=full_error)
        else:
            self.finish()

    # ----------------------------------
    # 公开 API
    # ----------------------------------

    def benchmark(self, name: str, config: Optional[Mapping[str, Any]] = None) -> Benchmark:
        """声明一个评测项。

        评测项名称即分数记录的 dataset 维度，评测项的局部配置并入
        ``config.benchmarks`` 并立即写入配置组件。

        :param name: 评测项名称，1..512 字符。
        :param config: 评测项的局部配置，例如 ``{"num_fewshot": 5}``。
        :return: 评测项对象，用于上报分数。
        :raises ValueError: 名称非法、配置无法序列化，或同名评测项已用不同配置声明。
        """
        if fork.is_forked(self._init_pid):
            raise RuntimeError(
                "SwanLab Eval does not support fork yet. Use `multiprocessing.set_start_method('spawn')` "
                "or call `swanlab.eval()` in the child process."
            )
        if not self.alive:
            raise RuntimeError("`run.benchmark()` requires an active Eval, call `swanlab.eval()` first.")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Benchmark name must be a non-empty string.")
        if len(name) > MAX_BENCHMARK_NAME_LENGTH:
            raise ValueError(f"Benchmark name exceeds the maximum length of {MAX_BENCHMARK_NAME_LENGTH} characters.")
        this_config = self._validate_benchmark_config(name, config)
        with self._api_lock:
            if not self.alive:
                raise RuntimeError("`run.benchmark()` requires an active Eval, call `swanlab.eval()` first.")
            existing = self._benchmarks.get(name)
            if existing is not None:
                if existing.config != this_config:
                    raise ValueError(f"Benchmark '{name}' has already been declared with a different config.")
                return existing
            benchmark = Benchmark(self, name, this_config)
            self._benchmarks[name] = benchmark
            self.config.update({"benchmarks": {item_name: item.config for item_name, item in self._benchmarks.items()}})
            return benchmark

    @staticmethod
    def _validate_benchmark_config(name: str, config: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
        if config is None:
            return {}
        if not isinstance(config, Mapping):
            raise ValueError(f"Benchmark config for '{name}' must be a dict, but got {type(config).__name__}.")
        try:
            return parse(dict(config))
        except TypeError as e:
            raise ValueError(f"Benchmark config for '{name}' is not serializable: {e}") from None

    def finish(self, state: FinishType = "success", error: Optional[str] = None) -> None:
        """Finish the current eval and wait for all scores to be flushed.

        :param state: Terminal state of the eval. Defaults to ``"success"``.
        :param error: Optional error message, required when ``state`` is ``"crashed"``.
        """
        if not self.alive:
            console.warning("SwanLab Eval has already finished or has not started.")
            return
        state = cast(FinishType, state.lower())
        if state not in _FINISH_STATES:
            console.error(f"Invalid state: {state}, allowed values are {sorted(_FINISH_STATES)}")
            return
        if state == "crashed" and error is None:
            console.warning("Crashed reason has been set to 'unknown' due to missing error message.")
            error = "unknown"
        with self._api_lock:
            if not self.alive:
                return
            # 1. 运行结束，停止组件，config 的缓冲记录在此刷入 core
            self._state = "finished"
            self._components.stop()
            # 2. 停止 Core：结束记录写入本地，第二阶段的确认只等待上传排空
            ts = Timestamp()
            ts.GetCurrentTime()
            finish_request = DeliverEvalFinishRequest(
                finish_record=EvalFinishRecord(
                    state=adapter.eval_state[state],
                    error=error,
                    finished_at=ts,
                )
            )
            finish_resp = self._core.deliver_eval_finish(finish_request)
            if not finish_resp.success:
                console.error(finish_resp.message)
            confirm_resp = run_with_progress(
                stats_fn=self._core.get_operation_stats,
                blocking_fn=self._core.confirm_eval_finish,
            )
            if not confirm_resp.success:
                console.error(confirm_resp.message)
            # 3. 清理副作用
            atexit.unregister(self._handle_atexit)
            sys.excepthook = self._sys_origin_excepthook
            if self._is_main_thread:
                signal.signal(signal.SIGINT, self._original_sigint_handler)
            clear_eval_run()
            if self.mode == "online" and client.exists():
                console.debug("Reset online client...")
                client.reset()
            console.debug("SwanLab Eval has finished.")
            console.reset()


_current_eval_run: Optional[EvalRun] = None


def has_eval_run() -> bool:
    """检查当前是否存在活跃的评估。

    :return: 是否存在活跃的 EvalRun。
    """
    return _current_eval_run is not None and _current_eval_run.alive


def get_eval_run() -> EvalRun:
    """获取当前活跃的评估。

    :return: 活跃的 EvalRun 实例。
    :raises RuntimeError: 当前没有活跃的评估。
    """
    if _current_eval_run is None:
        raise RuntimeError("No active Eval. Call swanlab.eval() first.")
    return _current_eval_run


def set_eval_run(run: EvalRun) -> None:
    global _current_eval_run
    _current_eval_run = run


def clear_eval_run() -> None:
    global _current_eval_run
    _current_eval_run = None
