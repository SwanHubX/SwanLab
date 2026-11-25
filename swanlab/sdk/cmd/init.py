"""
@author: cunyue
@file: init.py
@time: 2026/3/6 21:47
@description: SwanLab SDK 初始化方法
我们在设计上将init前后作为分界线，在init之前出现的错误（如配置校验失败、或 `on_online_error="raise"` 下的登录失败）被视为critical错误，直接抛出异常；
在online模式下默认会按 `on_online_error` 策略进行后台恢复或降级，在init之后的swanlab内部错误被视为non-critical错误，会尝试继续运行，但会记录错误日志
init函数执行时被视为init之前
"""

import json
import os
import socket
import time
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Any, Callable, Dict, List, NoReturn, Optional, Tuple, Union

import requests
import yaml
from google.protobuf.message import DecodeError
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.exceptions import AuthenticationError, DataStoreError
from swanlab.proto.swanlab.grpc.core.v1.core_pb2 import (
    CORE_ERROR_KIND_AUTH,
    CORE_ERROR_KIND_TRANSIENT,
    DeliverRunStartRequest,
    DeliverRunStartResponse,
)
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import StartRecord
from swanlab.sdk.cmd.guard import with_cmd_lock
from swanlab.sdk.internal.context import (
    RunConfig,
    RunContext,
    use_context,
)
from swanlab.sdk.internal.context.components import CoreBackend, create_core, create_probe, resolve_core_backend
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.store import DataStoreReader
from swanlab.sdk.internal.pkg import adapter, console, fork, fs, helper, safe
from swanlab.sdk.internal.pkg.client import is_login_auth_error, is_recoverable_login_error, is_transient_login_error
from swanlab.sdk.typings.cmd import ConfigLike
from swanlab.sdk.typings.context import CallbacksType
from swanlab.utils.experiment import generate_color, generate_id, generate_name

from ..internal.run import Run, get_run, has_run
from ..internal.settings import Settings, create_settings
from ..typings.run import ModeType, OnlineErrorPolicy, ParallelType, ResumeType
from . import utils
from .login import login_cli, login_raw

__all__ = [
    "init",
    "compatible_kwargs",
    "load_config",
    "prompt_init_mode",
    "set_nested_value",
    "ensure_run_dir",
]


def set_nested_value(d: dict, key: str, value: Any):
    """
    根据点分隔的键路径设置嵌套字典中的值
    例如: set_nested_value(d, "a.b.c", 1) 会设置 d["a"]["b"]["c"] = 1
    注意: 如果 value 为 None，则不会设置该值
    :param d: 要操作的字典
    :param key: 键路径，使用点分隔
    :param value: 要设置的值，如果为 None 则不设置
    """
    if value is None:
        return
    keys = key.split(".")
    current_dict = d
    for k in keys[:-1]:
        if k not in current_dict:
            current_dict[k] = {}
        current_dict = current_dict[k]
    current_dict[keys[-1]] = value


def compatible_kwargs(model_dict: dict, **kwargs) -> dict:
    """
    由于不同库的参数名不同，并且照顾到用户习惯，我们需要针对一些参数进行兼容性处理。
    将一些额外的参数合并到 model_dict 中
    """
    # experiment_name --> name
    set_nested_value(model_dict, "experiment.name", kwargs.pop("experiment_name", None))
    # notes --> description
    set_nested_value(model_dict, "experiment.description", kwargs.pop("notes", None))
    # logdir --> log_dir (backward compatibility)
    set_nested_value(model_dict, "log_dir", kwargs.pop("logdir", None))
    return model_dict


@with_cmd_lock
def init(
    *,
    reinit: Optional[bool] = None,
    log_dir: Optional[str] = None,
    mode: Optional[ModeType] = None,
    on_online_error: Optional[OnlineErrorPolicy] = None,
    workspace: Optional[str] = None,
    project: Optional[str] = None,
    public: Optional[bool] = None,
    name: Optional[str] = None,
    color: Optional[str] = None,
    description: Optional[str] = None,
    job_type: Optional[str] = None,
    group: Optional[str] = None,
    tags: Optional[List[str]] = None,
    id: Optional[str] = None,
    resume: Optional[Union[ResumeType, bool]] = None,
    parallel: Optional[ParallelType] = None,
    config: Optional[ConfigLike] = None,
    settings: Optional[Settings] = None,
    callbacks: Optional[CallbacksType] = None,
    **kwargs,
) -> Run:
    """Initialize a new SwanLab run to track experiments.

    This function starts a new run for logging metrics, artifacts, and metadata.
    After calling this, use `swanlab.log()` to log data and `swanlab.finish()` to
    close the run. SwanLab automatically finishes runs at program exit.

    :param reinit: If True, finish the current run before starting a new one. Defaults to False.
    :param log_dir: Directory to store logs. Defaults to "./swanlog".
    :param mode: Run mode. Options: "online" (sync to cloud), "local" (local only),
        "offline" (save locally for later sync), "disabled" (no logging). Defaults to "online".
    :param on_online_error: Behavior when SwanLab online service fails to connect or authenticate
        in online mode. Transient connection/server failures and authentication errors (401/403)
        are handled according to this policy; configuration errors still raise. Options: "offline"
        (keep online mode and save locally while retrying the connection on transient outages, or
        fall back to offline mode on auth/unrecoverable errors), "disabled" (no logging), "raise".
        Defaults to "offline".
    :param workspace: Workspace or organization name. Defaults to current user.
    :param project: Project name. Defaults to current directory name.
    :param public: Make project publicly visible (online mode only). Defaults to False.
    :param name: Experiment name. Auto-generated if not provided.
    :param color: Experiment color for visualization. Auto-generated if not provided.
    :param description: Experiment description.
    :param job_type: Job type label (e.g., "train", "eval").
    :param group: Group name for organizing related experiments.
    :param tags: List of tags for categorizing experiments.
    :param id: Run ID for identifying this run (used for cloud resume in online mode and
        local run identity verification when reusing ``run.dir``).
    :param resume: Cloud resume behavior. Options: "must" (must resume), "allow" (resume if exists),
        "never" (always create new). Defaults to "never". Only governs cloud experiment reuse;
        to resume a local run directory after a crash, explicitly set ``run.dir`` or ``SWANLAB_RUN_DIR``
        (requires a working Go core).
    :param parallel: Parallel execution strategy. Options: "none" (default), "shared".
        ``"shared"`` implies ``resume="allow"``, a hostname+pid suffixed run directory, and
        (unless ``core.detach`` is set explicitly) one detached ``swanlab-core`` shared by all ranks.
    :param config: Experiment configuration dict or path to config file (JSON/YAML).
    :param settings: Custom Settings object for advanced configuration.
    :param callbacks: List of callback functions triggered on run events.

    :return: The initialized Run object.

    :raises RuntimeError: If a run is already active and reinit=False, or if an online
        startup fails and is not degraded (a non-transient error, or a transient error
        when ``on_online_error="raise"``).
    :raises AuthenticationError: If ``on_online_error="raise"`` and SwanLab online service
        rejects the API key or login fails (wrapping the original exception in ``__cause__``),
        or if initial login fails with a non-recoverable non-auth error. Under the default
        ``on_online_error`` policies, 401/403 authentication errors degrade instead of raising.

    Examples:

        Basic local run:

        >>> import swanlab
        >>> swanlab.init(mode="local", project="my_project")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Online run with configuration:

        >>> import swanlab
        >>> swanlab.login(api_key="your_key")
        >>> swanlab.init(
        ...     mode="online",
        ...     project="image_classification",
        ...     name="resnet50_experiment",
        ...     config={"lr": 0.001, "batch_size": 32}
        ... )
        >>> swanlab.log({"accuracy": 0.95})
        >>> swanlab.finish()

        Resume a previous run:

        >>> import swanlab
        >>> swanlab.init(
        ...     mode="online",
        ...     project="my_project",
        ...     id="previous_run_id",
        ...     resume="must"
        ... )
    """
    if reinit and has_run():
        run = get_run()
        run.finish()
    if has_run():
        raise RuntimeError(
            "`swanlab.init` requires an inactive Run. Please use `swanlab.finish()` or `swanlab.init(reinit=True)` first."
        )
    # 运行时配置
    run_settings = create_settings()
    # --------------------- 合并配置，检查格式，与业务无关 ----------------------------
    # 配置具有优先级，从低到高依次是：全局配置 --> 自定义配置 --> 传入的参数
    # 1. 合并自定义配置
    if settings:
        run_settings.merge_settings(settings)
    # 2. 基于传入的参数，合并当前配置，此步骤必须在合并全局配置之后，因为传入的参数的优先级是高于全局配置的
    args_dict = compatible_kwargs({}, **kwargs)
    for key, value in {
        "log_dir": log_dir,
        "mode": mode,
        "on_online_error": on_online_error,
        "project.name": project,
        "project.workspace": workspace,
        "project.public": public,
        "experiment.name": name,
        "experiment.color": color,
        "experiment.description": description,
        "experiment.job_type": job_type,
        "experiment.group": group,
        "experiment.tags": tags,
        "run.resume": resume,
        "run.id": id,
        "run.config": Path(config) if isinstance(config, (str, os.PathLike)) else None,
        "run.parallel": parallel,
    }.items():
        set_nested_value(args_dict, key, value)
    run_settings.merge_settings(args_dict)
    # 仅在合并后的配置仍没有项目名时，使用当前目录名兜底
    if run_settings.project.name is None:
        run_settings.merge_settings({"project": {"name": Path.cwd().name}})
    # ---------------------------------- 再次确认参数 ----------------------------------
    # 根据交互式引导确定最终的模式
    mode = prompt_init_mode(run_settings)
    if run_settings.core.skip_store and mode != "online":
        # 交互式引导允许用户从 online 降级到 offline；skip_store 仅 online 合法，
        # 此处显式关闭并告警，避免 merge_settings 抛出裸 ValidationError，
        # 同时保证离线数据正常落盘（否则与 “Results will be saved locally” 提示矛盾）。
        console.warning("core.skip_store is only supported in online mode; it has been disabled for this run.")
        run_settings.merge_settings({"mode": mode, "core": {"skip_store": False}})
    else:
        run_settings.merge_settings({"mode": mode})
    # 校验 run id 与 resume，仅在对两者存在性有要求的模式下校验
    if run_settings.mode == "online":
        if run_settings.run.resume == "must":
            assert run_settings.run.id is not None, "Run id must be provided when resume=must."
        elif run_settings.run.resume == "never":
            assert run_settings.run.id is None, "Run id should not be provided when resume=never."
    # ---------------------------------- 初始化 ----------------------------------
    # 开始初始化
    generate_ctx = _init
    if run_settings.mode == "online":
        generate_ctx = utils.with_loading_animation()(_init)
    ctx, path, fallback_error = generate_ctx(run_settings, callbacks)
    # 初始化run
    try:
        run = Run(ctx, path, fallback_error=fallback_error)
    except BaseException as exc:
        # Run.__init__ 内部已回滚自身副作用；这里兜底回收构造早期失败时已启动的 Core
        with safe.block(message="core abort failed during init"):
            ctx.core.abort(exc)
        raise
    try:
        # 发送webhook回调，在除了disabled模式外，都会触发；尚未绑定云端路径时不拼接假链接
        success = send_webhook(ctx, path=path)
        if not success:
            webhook_url = ctx.config.settings.integration.webhook.url
            console.warning(
                f"Failed to send webhook, maybe due to network issues or invalid webhook URL: {webhook_url}."
            )
        # 加载配置并合并到 run.config
        config_data = load_config(run_settings, config)
        if config_data:
            run.config.update(config_data)
        return run
    except BaseException as exc:
        # Run 构造完成后的失败：回滚已启动的 Core、组件与全局 hook
        run._abort_init(exc)
        raise


def load_config(run_settings: Settings, config: Optional[ConfigLike]) -> Dict[str, Any]:
    """
    优雅地加载配置：支持字典直接返回，或从 JSON/YAML 文件加载。
    """
    config = config or run_settings.run.config
    # 1. 如果 config 为 None，则返回空字典
    if config is None:
        return {}

    # 2. 如果 config 是字典，则直接返回
    if isinstance(config, dict):
        return config

    # 3. 处理路径类型（包括字符串字面量）
    if isinstance(config, (str, os.PathLike)):
        path = Path(config)

        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        # 根据后缀名选择加载器
        suffix = path.suffix.lower()
        try:
            with open(path, "r", encoding="utf-8") as f:
                if suffix == ".json":
                    return json.load(f)
                elif suffix in (".yaml", ".yml"):
                    # 使用 safe_load 保证安全
                    return yaml.safe_load(f)
                else:
                    # 如果没有后缀或者后缀未知，尝试先按 JSON 读，失败再按 YAML 读
                    content = f.read()
                    try:
                        return json.loads(content)
                    except json.JSONDecodeError:
                        return yaml.safe_load(content)
        except Exception as e:
            raise ValueError(f"Error parsing config file {path}: {e}")

    # 4. 如果 config 不是上述类型，直接报错
    raise ValueError(f"Invalid config type: {type(config).__name__}. Expected dict, str, or PathLike.")


def prompt_init_mode(settings: Settings) -> ModeType:
    """
    在 swanlab.init 阶段，针对 online 模式且未登录的用户进行交互式引导。

    规则：
    1. 只有在 settings.interactive 为 True 且 settings.mode 为 'online' 时触发 。
    2. 如果 client 已存在（已登录），直接跳过 。
    3. 提供三个选项：(1) 使用已有的Key (2) 注册 (3) 切换为 offline 模式。

    :param settings: 当前的 Settings 实例 。
    :return: 最终确定的 mode
    """
    # 如果不是云模式，或者已经登录，或者非交互环境，直接返回当前状态
    mode = settings.mode
    if mode != "online" or client.exists():
        return mode
    login_func = partial(login_cli, save=True, host=settings.api_host)
    if mode == "online":
        if settings.api_key is not None:
            # 不登录，交给后面处理，否则会出现闪烁动画，比较影响美感
            # login_func(api_key=settings.api_key)
            return "online"

        if not settings.interactive:
            raise RuntimeError(
                "Failed to initialize SwanLab in online mode: no API key was provided, "
                "and interactive prompts are disabled."
            )
        if not helper.is_interactive():
            raise RuntimeError(
                "Failed to initialize SwanLab in online mode: no TTY is available for interactive login."
            )

        console.info("Using SwanLab to track your experiments. To get started, choose one of the following options:")
        console.print(
            "(1) Use an existing API key.",
            "(2) Create a new SwanLab account.",
            "(3) Continue without visualization (Offline mode).",
            "Learn more in the documentation: https://docs.swanlab.cn",
            sep="\n",
        )
        while True:
            choice = input("Enter your choice [1/2/3]: ").strip()

            if choice == "1":
                console.info("Using an existing SwanLab API key.")
                return _login_or_fallback(login_func, settings)

            if choice == "2":
                console.info(f"Create a SwanLab account here: {settings.web_host}/login")
                return _login_or_fallback(login_func, settings)

            if choice == "3":
                console.info("Continuing in Offline mode. Results will be saved locally.")
                return "offline"

            console.warning("Invalid choice. Please enter 1, 2, or 3.")
    # 其他模式不登录
    return mode


def _login_or_fallback(login_func: Callable[[], bool], settings: Settings) -> ModeType:
    """
    交互式引导中的登录：连接失败或认证失败（401/403）时按 on_online_error 降级，而不是终止训练进程。
    配置错误与 on_online_error="raise" 仍然抛错。
    """
    try:
        login_func()
    except Exception as e:
        is_auth = is_login_auth_error(e)
        if (not is_transient_login_error(e) and not is_auth) or settings.on_online_error == "raise":
            raise
        target_mode = settings.on_online_error
        _warn_fallback(target_mode, auth_error=is_auth)
        return target_mode
    return "online"


@safe.decorator(message="Failed to send webhook")
def send_webhook(ctx: RunContext, path: Optional[str]) -> Tuple[bool, bool]:
    """
    发送 webhook 回调，仅在非 disabled 模式下触发。

    请求体结构 (JSON):
    {
      "value": "string",  // 即 SWANLAB_WEBHOOK_VALUE 的值
      "swanlab": {
        "version": "string",     // swanlab 版本号
        "mode": "online" | "local", // swanlab 运行模式
        "run_dir": "string",  // 日志存储路径
        "exp_url": "string"       // 云端实验路径
      }
    }

    :param ctx: 运行上下文
    :param path: 已绑定的云端实验路径；为空表示尚未接入，不输出未生效的云端链接
    :return: (是否发送，是否成功)
    """
    if ctx.config.settings.mode == "disabled":
        console.debug("Skipping webhook because mode is disabled.")
        return False, False
    settings = ctx.config.settings
    webhook = settings.integration.webhook
    webhook_url = webhook.url
    if not webhook_url:
        console.debug("Skipping webhook because SWANLAB_WEBHOOK is not set.")
        return False, False
    webhook_value = webhook.value
    webhook_timeout = webhook.timeout
    # path 为空说明仍在等待后台恢复，不输出未生效的云端链接
    exp_url = f"{settings.web_host}{helper.fmt_run_path(path)}" if path else None
    # 发送请求
    requests.post(
        webhook_url,
        timeout=webhook_timeout,
        json={
            "value": webhook_value,
            "swanlab": {
                "version": helper.get_swanlab_version(),
                "mode": ctx.config.settings.mode,
                "run_dir": str(ctx.run_dir),
                "exp_url": exp_url,
            },
        },
    )
    return True, True


def ensure_run_dir(
    log_dir: Path,
    run_id: str,
    max_retries: int = 10,
    retry_interval: float = 1.0,
    dir_max_length: int = 255,
    dir_name: Optional[str] = None,
    parallel: ParallelType = "none",
    reuse_existing: bool = False,
) -> Path:
    """
    原子化地创建一个不与已有目录冲突的 run_dir，或在显式允许时复用已有目录。

    通过 ``mkdir(exist_ok=False)`` 将"检查不存在 + 创建目录"合并为一个原子操作：
    创建成功即拥有该目录；若抛出 FileExistsError，则等待 retry_interval 秒后，用新时间戳重试，最多重试 max_retries 次。

    当 ``dir_name`` 不为 None 时，直接使用指定目录名，不重试；默认仍要求目录不存在，
    仅当 ``reuse_existing=True`` 时允许复用已有目录。

    :param log_dir: 日志根目录（必须已存在）
    :param run_id: 当前运行的唯一标识符
    :param max_retries: 最大重试次数
    :param retry_interval: 目录冲突时等待的秒数
    :param dir_max_length: run_dir 目录名最大长度
    :param dir_name: 用户指定的 run_dir 目录名，指定后仅创建一次不重试
    :param parallel: 并行策略
    :param reuse_existing: 为 True 且指定 ``dir_name`` 时允许复用已存在目录
    :return: 创建或复用的 run_dir 路径
    :raises RuntimeError: 超过最大重试次数仍无法创建唯一目录
    :raises FileExistsError: 指定 ``dir_name`` 且 ``reuse_existing=False`` 时目录已存在，或目标路径为非目录文件
    """
    if dir_name is not None:
        run_dir = log_dir / dir_name
        fs.safe_mkdir(run_dir, ensure_clean=not reuse_existing)
        return run_dir
    for i in range(max_retries):
        run_dir_name, truncated = _generate_run_dir_name(run_id, dir_max_length, parallel)
        run_dir = log_dir / run_dir_name
        try:
            fs.safe_mkdir(run_dir, ensure_clean=True)
            if truncated:
                console.warning(f"Run directory name length exceeds limit, truncated to '{run_dir_name}'.")
            return run_dir
        except FileExistsError:
            if i == max_retries - 1:
                break
            time.sleep(retry_interval)
            console.debug(f"Retrying to create a unique run directory, attempt {i + 1} of {max_retries}")
    raise RuntimeError(f"Failed to create a unique run directory after {max_retries} attempts")


def _generate_run_dir_name(
    run_id: str, max_length: int, parallel: ParallelType = "none", hostname_max_len: int = 64
) -> Tuple[str, bool]:
    """
    生成 run_dir 的目录名，格式为 ``run-{timestamp}-{run_id}``，超长时截断。
    如果为共享模式或者分布式训练，增加hostname和pid两个参数，格式为 ``run-{timestamp}-{run_id}-{hostname}-{pid}``，
    语义为 "某一时刻、某个运行id、某个主机、某个进程开启的运行"。

    :param run_id: 当前运行的唯一标识符
    :param max_length: 目录名最大字节长度
    :param parallel: 并行策略
    :param hostname_max_len: 主机名最大字节长度，仅在共享模式或者分布式时使用，hostname_max_len 有可能会被进一步缩减以适应 max_length

    :return: run_dir 目录名（非完整路径）和是否截断
    """
    # 1. 先准备前后缀，因为只有run_id和hostname有可能被截断
    prefix = "-".join(["run", datetime.now().strftime("%Y%m%d_%H%M%S")]) + "-"
    suffix: Optional[str] = None
    hostname: Optional[str] = None
    # 并行模式添加hostname和pid到目录名
    if parallel != "none":
        with safe.block(level="debug", message="Failed to get hostname for run directory name, using fallback."):
            hostname = socket.gethostname() or "unknown_hostname"
        if not hostname:
            # 虽然console.warning有可能会造成大量重复日志刷屏，但这种情况确实很罕见，且通常意味着系统环境有问题，
            # 用户应该被告知这个潜在风险，所以我们选择在这里打warning而不是silent fail
            console.warning(
                "Failed to get hostname for run directory name, swanlab will use a fallback name.",
                "This may cause issues when running multiple processes on the same machine.",
                "Consider setting a custom run directory name (SWANLAB_RUN_DIR).",
            )
            hostname = "unknown_hostname"
        suffix = f"-{fork.current_pid()}"
    # 2. 计算中间部分的可用长度（算上"-"分隔符），并确定 hostname 的最大长度
    reserved_length = len(prefix.encode("utf-8")) + len((suffix or "").encode("utf-8"))
    available_length = max_length - reserved_length
    if available_length <= 0:
        # 由于 settings 中 dir_max_length 的约束，理论上不会出现 available_length <= 0 的情况，但这里做一个兜底，避免出现负数导致后续逻辑混乱
        raise ValueError(
            f"Failed to generate run directory name: max_length ({max_length}) is too small to accommodate prefix and suffix. "
            f"Consider increasing `dir_max_length` (SWANLAB_RUN_DIR_MAX_LENGTH) or specifying a custom `run.dir` (SWANLAB_RUN_DIR)."
        )
    # run id 最小会被截断为 "abc...yz" 样式，长度为 7，并且 run id 和 hostname 之间存在一个 "-" 分隔符，所以至少需要 8 个字符的长度才能保证 run id 的完整性
    if hostname is not None:
        hostname_max_len = min(hostname_max_len, available_length - 7 - 1)
    # 3. 开始截断和格式化中间部分的长度，hostname 和 run_id 共享 available_length 的长度，优先保证 hostname 的完整性
    hostname_truncated = False
    if hostname is not None:
        hostname = fs.safe_fmt(hostname, fallback="unknown_host")
        try:
            hostname, hostname_truncated = fs.safe_truncate(hostname, hostname_max_len)
        except ValueError:
            raise ValueError(
                f"Failed to generate run directory name: max_length ({max_length}) is too small to truncate hostname. "
                f"Consider increasing `dir_max_length` (SWANLAB_RUN_DIR_MAX_LENGTH) or specifying a custom `run.dir` (SWANLAB_RUN_DIR)."
            ) from None
        available_length -= len(hostname.encode("utf-8")) + 1  # 减去 hostname 的字节长度和一个 "-" 分隔符
    run_id = fs.safe_fmt(run_id, fallback="unknown_run")
    try:
        run_id, run_id_truncated = fs.safe_truncate(run_id, available_length)
    except ValueError:
        raise ValueError(
            f"Failed to generate run directory name: max_length ({max_length}) is too small to truncate run_id. "
            f"Consider increasing `dir_max_length` (SWANLAB_RUN_DIR_MAX_LENGTH) or specifying a custom `run.dir` (SWANLAB_RUN_DIR)."
        ) from None
    # 4. 组合最终的目录名
    dir_name = prefix + run_id
    if hostname is not None:
        dir_name += f"-{hostname}"
    if suffix is not None:
        dir_name += suffix
    return dir_name, hostname_truncated or run_id_truncated


def _validate_existing_run_dir(run_dir: Path, run_id: str, backend: CoreBackend) -> None:
    """校验待复用的 run_dir 内既有 run-*.swanlab 的首条 StartRecord 是否与当前 run_id 一致。"""
    run_files = sorted(run_dir.glob("run-*.swanlab"))
    if not run_files:
        return
    if len(run_files) > 1:
        raise RuntimeError(f"Multiple run-*.swanlab files found in '{run_dir}'; cannot determine which run to resume.")
    # 仅 Go core 支持从已有 run 文件续传；内置实现无法追加写入
    if backend.impl == "python":
        reason = f" ({backend.fallback_reason})" if backend.fallback_reason is not None else ""
        raise RuntimeError(
            f"Run directory '{run_dir}' already contains run file '{run_files[0].name}', "
            f"but local resume requires the Go core{reason}; choose a new run directory."
        )
    run_file = run_files[0]
    reader = DataStoreReader()
    opened = False
    try:
        reader.open(run_file)
        opened = True
        payload = reader.scan()
    except (AssertionError, DataStoreError, OSError) as exc:
        raise RuntimeError(f"Failed to read existing run file '{run_file}': {exc}") from exc
    finally:
        if opened:
            reader.close()

    if payload is None:
        return
    record = Record()
    try:
        record.ParseFromString(payload)
    except DecodeError as exc:
        raise RuntimeError(f"Failed to parse start record from '{run_file}': {exc}") from exc
    if not record.HasField("start"):
        raise RuntimeError(f"Existing run file '{run_file}' is missing a valid StartRecord as its first record.")
    existing_id = record.start.id
    if existing_id != run_id:
        raise RuntimeError(
            f"Run directory '{run_dir}' already belongs to another run "
            f"(existing id={existing_id!r}, current id={run_id!r})."
        )


def _init(
    run_settings: Settings, callbacks: Optional[CallbacksType]
) -> Tuple[RunContext, Optional[str], Optional[str]]:
    """
    初始化运行时配置，在这之前，所有引导式交互都已经完成
    上下文生命周期通过 `Run` 管理，而非全局 `ContextVar`
    """
    mode = run_settings.mode
    skip_store = run_settings.core.skip_store
    # 实验路径，/:username/:project_name/:slug(run_id)，与open api命名一致
    path = None
    fallback_error: Optional[str] = None
    # 1. 生成run_id
    if not run_settings.run.id:
        run_settings.merge_settings({"run": {"id": generate_id()}})
    run_id = run_settings.run.id
    assert run_id, "Run id is not provided."
    should_mkdirs = mode != "disabled" and not skip_store
    backend = resolve_core_backend(mode)
    # 2. 创建或复用运行目录
    if should_mkdirs:
        explicit_dir = run_settings.run.dir
        reusing_existing = explicit_dir is not None and (run_settings.log_dir / explicit_dir).is_dir()
        if reusing_existing:
            assert explicit_dir is not None
            _validate_existing_run_dir(run_settings.log_dir / explicit_dir, run_id, backend)
        # 安全创建目录，并写入 .gitignore（如果目录为空）
        helper.mkdir_and_append_gitignore(run_settings.log_dir)
        # 未显式指定 run.dir 时仍要求新建；仅显式指定且已存在时复用（不清理已有文件）
        run_dir = ensure_run_dir(
            run_settings.log_dir,
            run_id,
            max_retries=run_settings.run.dir_create_retries,
            dir_max_length=run_settings.run.dir_max_length,
            dir_name=explicit_dir,
            parallel=run_settings.run.parallel,
            reuse_existing=reusing_existing,
        )
        if reusing_existing:
            console.info(f"Reusing existing run directory '{run_dir}'; will continue the previous unfinished run.")
    else:
        # disabled / skip_store 下 run_dir 仅作逻辑路径，不创建目录
        if run_settings.run.dir:
            run_dir = run_settings.log_dir / run_settings.run.dir
        else:
            run_dir_name, _ = _generate_run_dir_name(run_id, run_settings.run.dir_max_length, run_settings.run.parallel)
            run_dir = run_settings.log_dir / run_dir_name
    # 3. 创建一个临时的上下文，避免出现任何问题导致上下文残留
    with use_context(
        RunContext(config=RunConfig(settings=run_settings, run_dir=run_dir), callbacks=callbacks, core_backend=backend)
    ) as ctx:
        assert run_settings.project.name, "Project name is required."
        # 在启动任何一次 Core 前固定恢复资格，避免试错过程中产生的文件影响判断
        allow_recover_online = _can_recover_online(backend, run_settings, run_dir)
        recover_online = False
        # 3.1. 前置操作
        # online 模式前置：确保 client 已就绪，连接或认证失败时按策略恢复或降级
        if mode == "online":
            fallback_mode, recover_online, fallback_error = _ensure_online_client_or_fallback(
                ctx,
                run_settings,
                run_id,
                run_dir,
                allow_recover_online=allow_recover_online,
            )
            if fallback_mode is not None:
                mode = fallback_mode
                should_mkdirs = mode == "offline"
            elif recover_online:
                _prepare_recover_online(run_settings, run_id, run_dir)
                should_mkdirs = True
        # local 模式前置：注入 LocalCallback 到 callbacks 中
        elif mode == "local":
            from swanlab.deprecated.local import LocalCallbacker

            ctx.callbacker.merge_callbacks(LocalCallbacker())
        # 3.2. 非云端模式，确定 name/color，合并到 settings
        if mode != "online":
            _ensure_offline_identity(run_settings, run_id)
        # 3.3. 统一调用 deliver_run_start（core 内部按 mode 分发：online 走网络，其余本地处理）
        start_request = _build_start_request(run_settings, run_id, run_dir, recover_online=recover_online)
        try:
            resp = ctx.core.deliver_run_start(start_request)
            if not resp.success:
                is_auth = resp.error_kind == CORE_ERROR_KIND_AUTH
                if (
                    mode != "online"
                    or (resp.error_kind != CORE_ERROR_KIND_TRANSIENT and not is_auth)
                    or run_settings.on_online_error == "raise"
                ):
                    _raise_start_failure(resp)
                if not is_auth and not recover_online and allow_recover_online:
                    # 具备恢复条件时保持 online 模式，带着恢复标记重拉 Go Core
                    recover_online = True
                    resp = _restart_with_recover_online(ctx, run_settings, run_id, run_dir, backend)
                    should_mkdirs = True
                else:
                    # 认证失败或无法后台恢复时，按策略切换到离线或禁用模式
                    mode = run_settings.on_online_error
                    if is_auth:
                        fallback_error = resp.message or "authentication failed"
                    resp = _restart_with_fallback(ctx, run_settings, run_id, run_dir, mode, auth_error=is_auth)
                    should_mkdirs = mode == "offline"
            path = resp.path or path
            if mode == "online" and not recover_online:
                assert path, "Initialization path failed when mode=online"
            ctx.next_step(user_step=resp.global_step)
            ctx.global_system_step = resp.global_system_step
            # 暂时不允许resume时、旧实验进行硬件信息采集、终端代理
            resume_limit_kwargs = {}
            if not resp.new_experiment:
                resume_limit_kwargs["probe.hardware"] = False
                resume_limit_kwargs["probe.requirements"] = False
                resume_limit_kwargs["probe.conda"] = False
                resume_limit_kwargs["probe.git"] = False
                resume_limit_kwargs["probe.swanlab"] = False
                resume_limit_kwargs["probe.monitor"] = False
                resume_limit_kwargs["console.proxy_type"] = "none"
                console.info(
                    "Hardware information collection, monitor, and terminal proxy have been disabled in resume mode."
                )
            # 3.4. 从 core 响应同步配置（online 模式会覆盖为服务端分配的值）
            # online 响应的 name 由服务端填写；offline / local / disabled 不填 name，此时回退到启动记录里的名字，避免被 strip_empty_str 过滤掉。
            experiment_name = resp.name or resp.run.name
            sync_args = {}
            merge_dict = helper.strip_none(
                {
                    "project.workspace": resp.run.workspace,
                    "project.name": resp.run.project,
                    "experiment.name": experiment_name,
                    "experiment.color": resp.run.color,
                    **resume_limit_kwargs,
                },
                strip_empty_str=True,
            )
            for key, value in merge_dict.items():
                set_nested_value(sync_args, key, value)
            run_settings.merge_settings(sync_args)
        except BaseException as exc:
            with safe.block(message="core abort failed during init"):
                ctx.core.abort(exc)
            raise
    # 4. 创建数据目录
    try:
        if should_mkdirs:
            fs.safe_mkdirs(ctx.media_dir, ctx.files_dir, ctx.debug_dir)
        return ctx, path, fallback_error
    except BaseException as exc:
        with safe.block(message="core abort failed during init"):
            ctx.core.abort(exc)
        raise


def _build_start_request(
    run_settings: Settings,
    run_id: str,
    run_dir: Path,
    *,
    recover_online: bool = False,
) -> DeliverRunStartRequest:
    """构建启动请求。降级或恢复重试时重新构建，确保使用最新的 settings。"""
    ts = Timestamp()
    ts.GetCurrentTime()
    start_record = StartRecord(
        project=run_settings.project.name,
        workspace=run_settings.project.workspace,
        public=run_settings.project.public,
        name=run_settings.experiment.name,
        color=run_settings.experiment.color,
        description=run_settings.experiment.description,
        job_type=run_settings.experiment.job_type,
        group=run_settings.experiment.group,
        tags=run_settings.experiment.tags,
        id=run_id,
        resume=adapter.resume[run_settings.run.resume],
        started_at=ts,
    )
    return DeliverRunStartRequest(
        start_record=start_record,
        core_settings=run_settings.to_core_proto(run_id=run_id, run_dir=run_dir, recover_online=recover_online),
    )


def _has_local_run_file(run_dir: Path) -> bool:
    """判断运行目录下是否已存在非空的 run-*.swanlab 本地记录文件。"""
    if not run_dir.is_dir():
        return False
    return any(p.is_file() and p.stat().st_size > 0 for p in run_dir.glob("run-*.swanlab"))


def _can_recover_online(backend: CoreBackend, run_settings: Settings, run_dir: Path) -> bool:
    """判断是否满足保持 online 模式并在后台持续恢复云端绑定的条件。"""
    return (
        backend.impl == "go"
        and run_settings.mode == "online"
        and run_settings.on_online_error == "offline"
        and bool(run_settings.api_key)
        and (run_settings.run.resume == "never" or _has_local_run_file(run_dir))
    )


def _prepare_recover_online(run_settings: Settings, run_id: str, run_dir: Path) -> None:
    """在后台恢复启动前关闭 skip_store、建好本地目录，并在全新运行时补齐名称与颜色。"""
    skip_store_disabled = run_settings.core.skip_store
    if skip_store_disabled:
        run_settings.merge_settings({"core": {"skip_store": False}})
    helper.mkdir_and_append_gitignore(run_settings.log_dir)
    fs.safe_mkdir(run_dir)
    # 复用已有本地文件时保留原有实验身份，只在新目录下生成默认名称和颜色
    if not _has_local_run_file(run_dir):
        _ensure_offline_identity(run_settings, run_id)
    _reset_client_after_fallback()
    suffix = " (core.skip_store has been disabled)" if skip_store_disabled else ""
    console.warning(
        "🚨 SwanLab online service is currently unavailable; records will be saved locally,"
        f" and connection will be retried in background{suffix}."
    )


def _deliver_start_or_raise(ctx: RunContext, request: DeliverRunStartRequest) -> DeliverRunStartResponse:
    """投递启动请求；失败时转换为用户可见异常。"""
    resp = ctx.core.deliver_run_start(request)
    if not resp.success:
        _raise_start_failure(resp)
    return resp


def _restart_with_recover_online(
    ctx: RunContext,
    run_settings: Settings,
    run_id: str,
    run_dir: Path,
    backend: CoreBackend,
) -> DeliverRunStartResponse:
    """在线启动遇临时故障且可恢复时，携带恢复开关重启 Core。"""
    _prepare_recover_online(run_settings, run_id, run_dir)
    with safe.block(level="debug", message="core abort failed before recover restart"):
        ctx.core.abort()
    ctx.core = create_core(ctx, backend=backend)
    ctx.probe = create_probe(core=ctx.core, disabled=False)
    return _deliver_start_or_raise(ctx, _build_start_request(run_settings, run_id, run_dir, recover_online=True))


def _switch_run_to_fallback(
    ctx: RunContext,
    run_settings: Settings,
    run_id: str,
    run_dir: Path,
    target_mode: ModeType,
    *,
    auth_error: bool = False,
) -> None:
    """切换降级模式并重建 Core/Probe，供登录降级与启动降级复用。"""
    _switch_run_mode(run_settings, target_mode, auth_error=auth_error)
    _ensure_offline_identity(run_settings, run_id)
    # 已知副作用：online 启动前已经创建了 log_dir 和空的 run_dir。
    # offline 降级会复用 run_dir；disabled 降级不会使用这两个目录，因此会遗留 log_dir（仅含 .gitignore）和空的 run_dir。
    # 清理需要判断 log_dir 是否由本次 init 创建，避免误删用户共享目录，收益不足以抵消复杂度，故接受该缺陷。
    if target_mode == "offline":
        # 原本 skip_store 时 run_dir 尚未创建，降级后必须补齐
        helper.mkdir_and_append_gitignore(run_settings.log_dir)
        fs.safe_mkdir(run_dir)
    # 先中止已启动的在线 Core（例如 Go Core 子进程），避免切换为本地 Core 时遗留孤儿进程
    with safe.block(level="debug", message="core abort failed before fallback switch"):
        ctx.core.abort()
    ctx.core = create_core(ctx)
    ctx.probe = create_probe(core=ctx.core, disabled=target_mode == "disabled")


def _restart_with_fallback(
    ctx: RunContext,
    run_settings: Settings,
    run_id: str,
    run_dir: Path,
    target_mode: ModeType,
    *,
    auth_error: bool = False,
) -> DeliverRunStartResponse:
    """在线启动失败后，按策略切换模式并在同一个 run_dir 上重新启动 Core。"""
    _switch_run_to_fallback(ctx, run_settings, run_id, run_dir, target_mode, auth_error=auth_error)
    return _deliver_start_or_raise(ctx, _build_start_request(run_settings, run_id, run_dir))


def _ensure_offline_identity(run_settings: Settings, run_id: str) -> None:
    """非 online 模式补齐实验名称与颜色。"""
    args_dict = {}
    for key, value in {
        "experiment.name": run_settings.experiment.name or generate_name("beauty"),
        "experiment.color": run_settings.experiment.color or generate_color("beauty"),
        "run.id": run_id,
    }.items():
        set_nested_value(args_dict, key, value)
    run_settings.merge_settings(args_dict)


def _warn_fallback(
    target_mode: ModeType,
    *,
    resume_must: bool = False,
    skip_store_disabled: bool = False,
    auth_error: bool = False,
) -> None:
    """确定降级策略后给出一次简短告警。上游错误已在其发生处打印。"""
    if target_mode == "offline":
        consequence = "data will be saved locally and can be synced later"
    else:
        consequence = "no data will be recorded"
    details = []
    if resume_must:
        details.append('resume="must" was not honored')
    if skip_store_disabled:
        details.append("core.skip_store has been disabled")
    suffix = f" ({'; '.join(details)})" if details else ""
    reason = (
        "Authentication with SwanLab online service failed" if auth_error else "SwanLab online service is unavailable"
    )
    console.warning(f"🚨 {reason}; falling back to {target_mode} mode: {consequence}{suffix}.")


def _reset_client_after_fallback() -> None:
    """降级后清理在线 client，避免下次 init 复用可能失效的 sid。"""
    with safe.block(message="Failed to reset SwanLab client after fallback"):
        if client.exists():
            console.debug("Resetting SwanLab client after fallback...")
            client.reset()


def _switch_run_mode(run_settings: Settings, target_mode: ModeType, *, auth_error: bool = False) -> None:
    """把本次 run 切换到目标模式，并完成告警与 client 清理。"""
    skip_store_disabled = run_settings.core.skip_store
    resume_must = run_settings.run.resume == "must"
    run_settings.merge_settings({"mode": target_mode, "core": {"skip_store": False}})
    _reset_client_after_fallback()
    _warn_fallback(
        target_mode,
        resume_must=resume_must,
        skip_store_disabled=skip_store_disabled,
        auth_error=auth_error,
    )


def _raise_start_failure(resp: DeliverRunStartResponse) -> NoReturn:
    """把 Core 返回的错误分类还原为用户可见异常。"""
    message = resp.message or "Failed to start run"
    if resp.error_kind == CORE_ERROR_KIND_AUTH:
        raise AuthenticationError(message)
    raise RuntimeError(message)


def _ensure_online_client_or_fallback(
    ctx: RunContext,
    run_settings: Settings,
    run_id: str,
    run_dir: Path,
    *,
    allow_recover_online: bool = False,
) -> Tuple[Optional[ModeType], bool, Optional[str]]:
    """
    确保 online client 就绪；连接失败或认证失败（401/403）时进入后台恢复或按 on_online_error 降级。
    只服务于 init 内部登录，显式 swanlab.login() 不受影响。

    :return: (降级后的模式, 是否进入后台恢复, 认证降级错误信息)
    """
    try:
        if not client.exists():
            assert run_settings.api_key, "API key is required."
            assert run_settings.api_host, "API host is required."
            login_raw(
                api_key=run_settings.api_key,
                host=run_settings.api_host,
                save=False,
                animation=False,
                print_welcome=False,
                quiet_transient=run_settings.on_online_error != "raise",
            )
        assert client.exists(), "No client found, please login first."
    except Exception as e:
        is_auth = is_login_auth_error(e)
        if (not is_transient_login_error(e) and not is_auth) or run_settings.on_online_error == "raise":
            raise
        if not is_auth and allow_recover_online and is_recoverable_login_error(e):
            return None, True, None
        target_mode = run_settings.on_online_error
        _switch_run_to_fallback(ctx, run_settings, run_id, run_dir, target_mode, auth_error=is_auth)
        fallback_msg = str(e.__cause__ or e) if is_auth else None
        return target_mode, False, fallback_msg
    return None, False, None
