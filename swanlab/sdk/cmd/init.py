"""
@author: cunyue
@file: init.py
@time: 2026/3/6 21:47
@description: SwanLab SDK 初始化方法
我们在设计上将init前后作为分界线，在init之前出现的错误（如登录失败）被视为critical错误，一旦报错直接退出（大多数情况下）
在init之后的swanlab内部错误被视为non-critical错误，会尝试继续运行，但会记录错误日志
init函数执行时被视为init之前
"""

import json
import os
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import requests
import yaml

from swanlab.sdk.cmd.guard import with_cmd_lock
from swanlab.sdk.internal.context import (
    RunConfig,
    RunContext,
    callbacker,
    use_context,
)
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.api.project import get_or_create_project, get_project
from swanlab.sdk.internal.pkg import console, fs, helper, safe
from swanlab.sdk.internal.protocol import Callback
from swanlab.utils import generate_color, generate_id, generate_name

from ..internal.core_python.api.experiment import create_or_resume_experiment
from ..internal.run import Run, get_run, has_run
from ..internal.settings import Settings
from ..internal.settings import settings as global_settings
from ..typings.run import ModeType, ResumeType
from . import utils
from .login import login_cli, login_raw

__all__ = ["init", "ConfigLike"]


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
    return model_dict


ConfigLike = Union[Dict[str, Any], str, os.PathLike]


@with_cmd_lock
def init(
    *,
    reinit: Optional[bool] = None,
    logdir: Optional[str] = None,
    mode: Optional[ModeType] = None,
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
    config: Optional[ConfigLike] = None,
    settings: Optional[Settings] = None,
    callbacks: Optional[List[Callback]] = None,
    **kwargs,
) -> Run:
    """Initialize a new SwanLab run to track experiments.

    This function starts a new run for logging metrics, artifacts, and metadata.
    After calling this, use `swanlab.log()` to log data and `swanlab.finish()` to
    close the run. SwanLab automatically finishes runs at program exit.

    :param reinit: If True, finish the current run before starting a new one. Defaults to False.

    :param logdir: Directory to store logs. Defaults to "./swanlog".

    :param mode: Run mode. Options: "cloud" (sync to cloud), "local" (local only),
        "offline" (save locally for later sync), "disabled" (no logging). Defaults to "cloud".

    :param workspace: Workspace or organization name. Defaults to current user.

    :param project: Project name. Defaults to current directory name.

    :param public: Make project publicly visible (cloud mode only). Defaults to False.

    :param name: Experiment name. Auto-generated if not provided.

    :param color: Experiment color for visualization. Auto-generated if not provided.

    :param description: Experiment description.

    :param job_type: Job type label (e.g., "train", "eval").

    :param group: Group name for organizing related experiments.

    :param tags: List of tags for categorizing experiments.

    :param id: Run ID for resuming a previous run (cloud mode only).

    :param resume: Resume behavior. Options: "must" (must resume), "allow" (resume if exists),
        "never" (always create new). Defaults to "never".

    :param config: Experiment configuration dict or path to config file (JSON/YAML).

    :param settings: Custom Settings object for advanced configuration.

    :param callbacks: List of callback functions triggered on run events.

    :return: The initialized Run object.

    :raises RuntimeError: If a run is already active and reinit=False.

    Examples:

        Basic local run:

        >>> import swanlab
        >>> swanlab.init(mode="local", project="my_project")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Cloud run with configuration:

        >>> import swanlab
        >>> swanlab.login(api_key="your_key")
        >>> swanlab.init(
        ...     mode="cloud",
        ...     project="image_classification",
        ...     name="resnet50_experiment",
        ...     config={"lr": 0.001, "batch_size": 32}
        ... )
        >>> swanlab.log({"accuracy": 0.95})
        >>> swanlab.finish()

        Resume a previous run:

        >>> import swanlab
        >>> swanlab.init(
        ...     mode="cloud",
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
    run_settings = Settings()
    # --------------------- 合并配置，检查格式，与业务无关 ----------------------------
    # 配置具有优先级，从低到高依次是：全局配置 --> 自定义配置 --> 传入的参数
    # 1. 合并全局配置
    run_settings.merge_settings(global_settings)
    # 2. 合并自定义配置
    if settings:
        run_settings.merge_settings(settings)
    # 3. 基于传入的参数，合并当前配置，此步骤必须在合并全局配置之后，因为传入的参数的优先级是高于全局配置的
    args_dict = compatible_kwargs({}, **kwargs)
    for key, value in {
        "logdir": logdir,
        "mode": mode,
        "project.name": project or Path.cwd().name,
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
    }.items():
        set_nested_value(args_dict, key, value)
    run_settings.merge_settings(args_dict)
    # ---------------------------------- 再次确认参数 ----------------------------------
    # 根据交互式引导确定最终的模式
    mode = prompt_init_mode(run_settings)
    run_settings.merge_settings({"mode": mode})
    # 校验 run id 与 resume，仅在对两者存在性有要求的模式下校验
    if run_settings.mode == "cloud":
        if run_settings.run.resume == "must":
            assert run_settings.run.id is not None, "Run id must be provided when resume=must."
        elif run_settings.run.resume == "never":
            assert run_settings.run.id is None, "Run id should not be provided when resume=never."
    # ---------------------------------- 初始化 ----------------------------------
    # 合并回调
    callbacker.merge_callbacks(callbacks or [])
    # 开始初始化
    ctx = _init(run_settings)
    # 初始化run
    run = Run(ctx)
    # 发送webhook回调，在除了disabled模式外，都会触发
    success = send_webhook(ctx)
    if not success:
        webhook_url = ctx.config.settings.integration.webhook.url
        console.warning(f"Failed to send webhook, maybe due to network issues or invalid webhook URL: {webhook_url}.")
    # 加载配置并合并到 run.config
    config_data = load_config(run_settings, config)
    if config_data:
        run.config.update(config_data)
    return run


def _init(run_settings: Settings) -> RunContext:
    """
    初始化运行时配置，在这之前，所有引导式交互都已经完成
    上下文生命周期通过 `Run` 管理，而非全局 `ContextVar`
    """
    mode = run_settings.mode
    # 生成run_id
    run_id = run_settings.run.id or generate_id()
    run_dir = run_settings.log_dir / ("run-" + datetime.now().strftime("%Y%m%d_%H%M%S") + "-" + run_id)
    # 创建一个临时的上下文，避免出现任何问题导致上下文残留
    with use_context(RunContext(config=RunConfig(settings=run_settings, run_dir=run_dir))) as ctx:
        assert run_settings.project.name, "Project name is required."
        # 根据模式进行特定处理
        if mode == "cloud":
            _init_cloud(ctx, run_id)
        elif mode == "local":
            _mkdirs(ctx)
            name = generate_name("beauty")
            color = generate_color("beauty")
            workspace = "local"
            # 合并配置
            args_dict = {}
            for key, value in {
                "experiment.name": name,
                "experiment.color": color,
                "run.id": run_id,
                "project.workspace": workspace,
            }.items():
                set_nested_value(args_dict, key, value)
            run_settings.merge_settings(args_dict)

            # TODO: 注册回调器
        elif mode == "offline":
            _mkdirs(ctx)
            name = generate_name("beauty")
            color = generate_color("beauty")
            workspace = "offline"
            # 合并配置
            args_dict = {}
            for key, value in {
                "experiment.name": name,
                "experiment.color": color,
                "project.workspace": workspace,
                "run.id": run_id,
            }.items():
                set_nested_value(args_dict, key, value)
            run_settings.merge_settings(args_dict)
        elif mode == "disabled":
            name = generate_name("beauty")
            color = generate_color("beauty")
            workspace = "disabled"
            # 合并配置
            args_dict = {}
            for key, value in {
                "experiment.name": name,
                "experiment.color": color,
                "project.workspace": workspace,
                "run.id": run_id,
            }.items():
                set_nested_value(args_dict, key, value)
            run_settings.merge_settings(args_dict)
            # 不注册回调器
            pass
        else:
            raise ValueError(f"Invalid mode for `swanlab.init`: {mode}")
    return ctx


@helper.with_loading_animation()
def _init_cloud(ctx: RunContext, run_id: str):
    """
    在云模式下初始化运行上下文。
    :param ctx: 运行上下文
    :param run_id: 当前运行的唯一标识符
    """
    run_settings = ctx.config.settings
    if not client.exists():
        assert run_settings.api_key, "API key is required."
        assert run_settings.api_host, "API host is required."
        login_raw(
            api_key=run_settings.api_key,
            host=run_settings.api_host,
            save=False,
            animation=False,
            wellcome_on_success=False,
        )
    assert run_settings.project.name, "Project name is required."
    assert client.exists(), "No client found, please login first."
    _mkdirs(ctx)
    # 获取当前项目，如果不存在则创建
    project = get_or_create_project(
        username=run_settings.project.workspace,
        name=run_settings.project.name,
        public=run_settings.project.public,
    )
    username, project = project["username"], project["name"]
    # 获取当前项目详细信息
    project_info = get_project(username=username, name=project)
    # 获取当前实验
    experiment = run_settings.experiment
    history_experiment_count = project_info["_count"]["experiments"]
    name = experiment.name or generate_name(history_experiment_count)
    color = experiment.color or generate_color(history_experiment_count)
    # 开启实验
    _ = create_or_resume_experiment(
        username,
        project,
        name=name,
        resume=run_settings.run.resume,
        run_id=run_id,
        color=color,
        description=experiment.description,
        job_type=experiment.job_type,
        group=experiment.group,
        tags=experiment.tags,
    )
    # TODO resume 时向后端获取数据或向本地获取数据

    # 最后同步一次配置
    args_dict = {}
    for key, value in {
        "project.workspace": username,
        "project.name": project,
        "experiment.name": name,
        "experiment.color": color,
        "run.id": run_id,
    }.items():
        set_nested_value(args_dict, key, value)
    run_settings.merge_settings(args_dict)


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


def _mkdirs(ctx: RunContext):
    """
    创建运行所需的目录
    :param ctx: 运行上下文
    """
    # 对于 logdir 而言，如果不存在则创建，如果为空则写入 .gitignore
    log_dir = ctx.config.settings.log_dir
    # 1. 安全创建目录（如果不存在）
    fs.safe_mkdir(log_dir)
    # 2. 写入 .gitignore（如果目录为空）
    utils.append_gitignore(log_dir)
    # 3. 创建别的目录
    fs.safe_mkdirs(ctx.run_dir, ctx.media_dir, ctx.files_dir, ctx.debug_dir)


def prompt_init_mode(settings: Settings) -> ModeType:
    """
    在 swanlab.init 阶段，针对 cloud 模式且未登录的用户进行交互式引导。

    规则：
    1. 只有在 settings.interactive 为 True 且 settings.mode 为 'cloud' 时触发 。
    2. 如果 client 已存在（已登录），直接跳过 。
    3. 提供三个选项：(1) 使用已有的Key (2) 注册 (3) 切换为 offline 模式。

    :param settings: 当前的 Settings 实例 。
    :return: (最终确定的 mode, 是否成功登录)
    """
    # 如果不是云模式，或者已经登录，或者非交互环境，直接返回当前状态
    mode = settings.mode
    if mode != "cloud" or client.exists() or not settings.interactive:
        return mode
    login_func = partial(login_cli, save=True, host=settings.api_host)
    if mode == "cloud":
        if settings.api_key is not None:
            # 不登录，交给后面处理，否则会出现闪烁动画，比较影响美感
            # login_func(api_key=settings.api_key)
            return "cloud"

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
                login_func(save=True)
                return "cloud"

            if choice == "2":
                console.info(f"Create a SwanLab account here:{settings.web_host}/login")
                login_func(save=True)
                return "cloud"

            if choice == "3":
                console.info("Continuing in Offline mode. Results will be saved locally.")
                return "offline"

            console.warning("Invalid choice. Please enter 1, 2, or 3.")
    # 其他模式不登录
    return mode


@safe.decorator(message="Failed to send webhook")
def send_webhook(ctx: RunContext) -> Tuple[bool, bool]:
    """
    发送 webhook 回调，仅在非 disabled 模式下触发。

    请求体结构 (JSON):
    {
      "value": "string",  // 即 SWANLAB_WEBHOOK_VALUE 的值
      "swanlab": {
        "version": "string",     // swanlab 版本号
        "mode": "cloud" | "local", // swanlab 运行模式
        "run_dir": "string",  // 日志存储路径
        "exp_url": "string"       // 云端实验路径
      }
    }

    :param ctx: 运行上下文
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
    # 获取实验url
    if settings.mode == "cloud":
        exp_url = f"{settings.web_host}/@{settings.project.workspace}/{settings.project.name}/runs/{settings.run.id}"
    else:
        exp_url = None
    # 发送请求
    requests.post(
        webhook_url,
        timeout=webhook_timeout,
        json={
            "value": webhook_value,
            "swanlab": {
                "version": helper.get_swanlab_version(),
                "mode": ctx.config.settings.mode,
                "run_dir": ctx.run_dir,
                "exp_url": exp_url,
            },
        },
    )
    return True, True
