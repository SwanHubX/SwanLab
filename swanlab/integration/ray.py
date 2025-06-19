"""
import random
from ray import tune
from swanlab.integration.ray import SwanLabLoggerCallback

def train_func(config):
    offset = random.random() / 5
    for epoch in range(2, config["epochs"]):
        acc = 1 - (2 + config["lr"]) ** -epoch - random.random() / epoch - offset
        loss = (2 + config["lr"]) ** -epoch + random.random() / epoch + offset
        tune.report({"acc": acc, "loss": loss})


tuner = tune.Tuner(
    train_func,
    param_space={
        "lr": tune.grid_search([0.001, 0.01, 0.1, 1.0]),
        "epochs": 10,
    },
    run_config=tune.RunConfig(
        callbacks=[SwanLabLoggerCallback(project="Ray_Project")],
    ),
)
results = tuner.fit()
"""

import enum
import os
import urllib
import warnings
from numbers import Number
from types import ModuleType
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pyarrow.fs

try:
    import ray
except ImportError:
    raise ImportError("Ray is not installed. Please install it using 'pip install ray'.")

from ray import logger
from ray.air.constants import TRAINING_ITERATION
from ray.air.util.node import _force_on_current_node
from ray.train._internal.session import get_session
from ray.train._internal.syncer import DEFAULT_SYNC_TIMEOUT
from ray.tune.experiment import Trial
from ray.tune.logger import LoggerCallback
from ray.tune.utils import flatten_dict
from ray.util import PublicAPI
from ray.util.queue import Queue

import swanlab


# 环境变量常量
SWANLAB_API_KEY = "SWANLAB_API_KEY"
SWANLAB_MODE = "SWANLAB_MODE"
SWANLAB_PROJ_NAME = "SWANLAB_PROJ_NAME"

@PublicAPI(stability="alpha")
def setup_swanlab(
    config: Optional[Dict] = None,
    api_key: Optional[str] = None,
    api_key_file: Optional[str] = None,
    swanlab_host: Optional[str] = None,
    rank_zero_only: bool = True,
    **kwargs,
):
    """
    设置SwanLab运行实例
    
    这是主要的初始化函数，用于在Ray训练会话中设置SwanLab日志记录。
    
    调用逻辑：
    1. 检查SwanLab是否已安装
    2. 获取当前训练会话信息（trial_id, trial_name, experiment_name）
    3. 调用_setup_swanlab进行实际初始化
    
    Args:
        config: 实验配置字典
        api_key: SwanLab API密钥
        api_key_file: 包含API密钥的文件路径
        rank_zero_only: 是否只在rank 0进程上初始化
        **kwargs: 传递给swanlab.init的其他参数
    
    Returns:
        SwanLab运行实例
    """
    # 默认值初始化
    default_trial_id = None
    default_trial_name = None
    default_experiment_name = None

    # 尝试获取当前训练会话信息
    # 如果不在训练会话中，session将为None
    session = get_session()
    if session and rank_zero_only and session.world_rank in (None, 0):
        return

    # 从会话中获取默认值
    if session:
        default_trial_id = session.trial_id
        default_trial_name = session.trial_name
        default_experiment_name = session.experiment_name

    # 构建SwanLab初始化参数
    # 优先级：kwargs > 默认值
    swanlab_init_kwargs = {
        "trial_id": kwargs.get("trial_id") or default_trial_id,
        "trial_name": kwargs.get("trial_name") or default_trial_name,
    }
    # 传入的kwargs覆盖默认kwargs
    swanlab_init_kwargs.update(kwargs)

    # 调用内部设置函数
    return _setup_swanlab(
        config=config, api_key=api_key, api_key_file=api_key_file, swanlab_host=swanlab_host, **swanlab_init_kwargs
    )


def _get_swanlab_project(project: Optional[str] = None) -> Optional[str]:
    """Get SwanLab project from environment variable or external hook if not passed
    as and argument."""
    if not project and os.environ.get(SWANLAB_PROJ_NAME):
        # Try to get project and group from environment variables if not
        # passed through SwanLabLoggerCallback.
        project = os.environ.get(SWANLAB_PROJ_NAME)
    return project


def _setup_swanlab(
    trial_id: str,
    trial_name: str,
    config: Optional[Dict] = None,
    api_key: Optional[str] = None,
    api_key_file: Optional[str] = None,
    swanlab_host: Optional[str] = None,
    _swanlab: Optional[ModuleType] = None,
    **kwargs,
):
    """
    内部SwanLab设置函数
    
    这是实际的SwanLab初始化函数，负责：
    1. 清理配置数据（移除不可序列化的对象）
    2. 设置API密钥
    3. 调用swanlab.init()创建运行实例
    
    调用逻辑：
    setup_swanlab -> _setup_swanlab -> swanlab.init()
    
    Args:
        trial_id: 试验ID
        trial_name: 试验名称
        config: 实验配置
        api_key: API密钥
        api_key_file: API密钥文件路径
        swanlab_host: SwanLab主机地址
        _swanlab: SwanLab模块（用于测试）
        **kwargs: 其他初始化参数
    
    Returns:
        SwanLab运行实例
    """
    _config = config.copy() if config else {}

    # 如果指定了密钥文件，展开用户路径
    # if api_key_file:
    #     api_key_file = os.path.expanduser(api_key_file)
    if api_key:
        swanlab.login(api_key=api_key, host=swanlab_host)

    project = _get_swanlab_project(kwargs.pop("project", None))

    # 移除不可序列化的项目
    _config = _clean_log(_config)

    # 构建SwanLab初始化参数
    swanlab_init_kwargs = dict(
        project=project,
        name=trial_name,
        reinit=True,  # 允许重新初始化
        config=_config,
    )

    # 更新配置（例如设置swanlab.init调用中的其他参数）
    swanlab_init_kwargs.update(**kwargs)

    _swanlab = _swanlab or swanlab

    # 初始化SwanLab运行实例
    run = _swanlab.init(**swanlab_init_kwargs)
    return run


def _is_allowed_type(obj):
    """
    检查对象类型是否允许记录到SwanLab
    
    只允许数字类型和大小为1的numpy数组（其元素为数字）
    
    Args:
        obj: 要检查的对象
    
    Returns:
        bool: 如果类型允许则返回True
    """
    if isinstance(obj, np.ndarray) and obj.size == 1:
        return isinstance(obj.item(), Number)
    return isinstance(obj, (Number))


def _clean_log(obj: Any):
    """
    清理日志对象，移除不可序列化的项目
    递归处理字典、列表、集合和元组，只保留允许的类型
    
    Args:
        obj: 要清理的对象
    
    Returns:
        清理后的对象
    """
    if isinstance(obj, dict):
        return {k: _clean_log(v) for k, v in obj.items()}
    elif isinstance(obj, (list, set)):
        return [_clean_log(v) for v in obj]
    elif isinstance(obj, tuple):
        return tuple(_clean_log(v) for v in obj)
    elif _is_allowed_type(obj):
        return obj
    
    fallback = str(obj)
    try:
        fallback = int(fallback)
        return fallback
    except ValueError:
        pass
    return fallback


class _QueueItem(enum.Enum):
    """
    队列项目类型枚举
    
    用于在SwanLab日志记录Actor中标识不同类型的消息
    """
    END = enum.auto()  # 结束信号
    RESULT = enum.auto()  # 结果数据
    CHECKPOINT = enum.auto()  # 检查点数据


class _swanlabLoggingActor:
    """
    SwanLab日志记录Actor类
    
    这是一个Ray Actor，在后台运行并负责：
    1. 初始化SwanLab运行实例
    2. 从队列接收日志数据
    3. 将数据发送到SwanLab服务器
    4. 处理配置更新
    
    调用逻辑：
    SwanLabLoggerCallback._start_logging_actor -> 创建Actor实例 -> Actor.run()
    SwanLabLoggerCallback.log_trial_result -> 发送数据到队列 -> Actor处理数据
    """
    
    def __init__(
        self,
        logdir: str,
        queue: Queue,
        exclude: List[str],
        to_config: List[str],
        *args,
        **kwargs,
    ):
        """
        初始化日志记录Actor
        
        Args:
            logdir: 日志目录
            queue: 用于接收日志数据的队列
            exclude: 要排除的结果键列表
            to_config: 要添加到配置的键列表
            *args, **kwargs: 传递给swanlab.init的参数
        """
        import swanlab

        self._swanlab = swanlab

        # 切换到日志目录
        os.chdir(logdir)
        self.queue = queue
        self._exclude = set(exclude)
        self._to_config = set(to_config)
        self.args = args
        self.kwargs = kwargs

        self._trial_name = self.kwargs.get("name", "unknown")
        self._logdir = logdir

    def run(self):
        """
        Actor的主运行循环
        
        调用逻辑：
        1. 初始化SwanLab运行实例
        2. 进入无限循环，从队列接收消息
        3. 根据消息类型处理数据
        4. 收到END信号时退出循环
        5. 调用swanlab.finish()完成运行
        """
        # 初始化SwanLab运行实例
        run = self._swanlab.init(*self.args, **self.kwargs)
        run.config.trial_log_path = self._logdir

        # 主处理循环
        while True:
            item_type, item_content = self.queue.get()
            if item_type == _QueueItem.END:
                break

            if item_type == _QueueItem.CHECKPOINT:
                # 处理检查点（当前未实现）
                # self._handle_checkpoint(item_content)
                continue

            # 处理结果数据
            assert item_type == _QueueItem.RESULT
            
            log, config_update = self._handle_result(item_content)
            try:
                # 更新配置并记录日志
                self._swanlab.config.update(config_update, allow_val_change=True)
                self._swanlab.log(log, step=log.get(TRAINING_ITERATION))
            except urllib.error.HTTPError as e:
                logger.warning("Failed to log result to swanlab: {}".format(str(e)))
                
        # 完成SwanLab运行
        self._swanlab.finish()
        
    def _handle_result(self, result: Dict) -> Tuple[Dict, Dict]:
        """
        处理结果数据
        
        将结果数据分离为日志数据和配置更新数据
        
        Args:
            result: 原始结果字典
        
        Returns:
            Tuple[Dict, Dict]: (日志数据, 配置更新数据)
        """
        config_update = result.get("config", {}).copy()
        log = {}
        # 展平结果字典
        flat_result = flatten_dict(result, delimiter="/")

        # 处理每个键值对
        for k, v in flat_result.items():
            # 跳过排除的键
            if any(k.startswith(item + "/") or k == item for item in self._exclude):
                continue
            # 添加到配置的键
            elif any(k.startswith(item + "/") or k == item for item in self._to_config):
                config_update[k] = v
            # 跳过不允许的类型
            elif not _is_allowed_type(v):
                continue
            # 添加到日志
            else:
                log[k] = v

        # 移除callbacks（不可序列化）
        config_update.pop("callbacks", None)
        return log, config_update


@PublicAPI(stability="alpha")
class SwanLabLoggerCallback(LoggerCallback):
    """
    SwanLab日志记录回调类
    
    这是Ray Tune的主要回调类，负责：
    1. 为每个试验创建和管理日志记录Actor
    2. 处理试验生命周期事件（开始、结果、结束）
    3. 协调多个试验的日志记录
    
    调用逻辑：
    Ray Tune -> SwanLabLoggerCallback.log_trial_start -> _start_logging_actor
    Ray Tune -> SwanLabLoggerCallback.log_trial_result -> 发送数据到Actor队列
    Ray Tune -> SwanLabLoggerCallback.log_trial_end -> 停止Actor并清理
    """
    
    # 不记录这些结果键
    _exclude_results = ["done", "should_checkpoint"]

    # 自动添加到配置的键
    AUTO_CONFIG_KEYS = [
        "trial_id",
        "experiment_tag",
        "node_ip",
        "experiment_id",
        "hostname",
        "pid",
        "date",
    ]

    _logger_actor_cls = _swanlabLoggingActor

    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        api_key_file: Optional[str] = None,
        api_key: Optional[str] = None,
        swanlab_host: Optional[str] = None,
        excludes: Optional[List[str]] = None,
        log_config: bool = False,
        upload_checkpoints: bool = False,
        save_checkpoints: bool = False,
        upload_timeout: int = DEFAULT_SYNC_TIMEOUT,
        **kwargs,
    ):
        """
        初始化SwanLab日志记录回调
        
        Args:
            project: SwanLab项目名称
            group: 实验组名称
            api_key_file: API密钥文件路径
            api_key: API密钥
            excludes: 要排除的结果键列表
            log_config: 是否在每个结果中记录配置
            upload_checkpoints: 是否上传检查点
            save_checkpoints: 已弃用，使用upload_checkpoints
            upload_timeout: 上传超时时间
            **kwargs: 其他参数
        """
        if not swanlab:
            raise RuntimeError(
                "swanlab was not found - please install with `pip install swanlab`"
            )

        # 处理弃用参数
        if save_checkpoints:
            warnings.warn(
                "`save_checkpoints` is deprecated. Use `upload_checkpoints` instead.",
                DeprecationWarning,
            )
            upload_checkpoints = save_checkpoints

        # 初始化属性
        self.project = project
        self.workspace = workspace
        self.api_key_path = api_key_file
        self.api_key = api_key
        self.swanlab_host = swanlab_host
        self.excludes = excludes or []
        self.log_config = log_config
        self.upload_checkpoints = upload_checkpoints
        self._upload_timeout = upload_timeout
        self.kwargs = kwargs

        # Actor相关属性
        self._remote_logger_class = None

        # 试验管理字典
        self._trial_logging_actors: Dict[
            "Trial", ray.actor.ActorHandle[_swanlabLoggingActor]
        ] = {}
        self._trial_logging_futures: Dict["Trial", ray.ObjectRef] = {}
        self._logging_future_to_trial: Dict[ray.ObjectRef, "Trial"] = {}
        self._trial_queues: Dict["Trial", Queue] = {}

    def setup(self, *args, **kwargs):
        """
        设置回调
        
        在实验开始时调用，进行初始化设置
        """
        # if self.api_key:
        #     swanlab.login(api_key=self.api_key, host=self.swanlab_host)

        self.project = _get_swanlab_project(self.project)
        if not self.project:
            raise ValueError(
                "Please pass the project name as argument or through "
                f"the {SWANLAB_PROJ_NAME} environment variable."
            )

    def log_trial_start(self, trial: "Trial"):
        """
        记录试验开始
        
        当新试验开始时调用，负责：
        1. 准备试验配置
        2. 创建日志记录Actor
        3. 启动后台日志记录
        
        调用逻辑：
        Ray Tune -> log_trial_start -> _start_logging_actor -> 创建Actor
        
        Args:
            trial: Ray Tune试验对象
        """
        config = trial.config.copy()

        # 移除callbacks（不可序列化）
        config.pop("callbacks", None)

        # 构建排除列表
        exclude_results = self._exclude_results.copy()
        exclude_results += self.excludes

        # 如果不记录配置，排除config键
        if not self.log_config:
            exclude_results += ["config"]

        # 获取试验ID和名称
        trial_id = trial.trial_id if trial else None
        trial_name = str(trial) if trial else None

        # SwanLab项目名称和工作空间
        swanlab_project = self.project
        swanlab_workspace = self.workspace

        # 清理配置（移除不可序列化的项目）
        config = _clean_log(config)
        config = {
            key: value for key, value in config.items() if key not in self.excludes
        }

        # 构建SwanLab初始化参数
        swanlab_init_kwargs = dict(
            project=swanlab_project,
            workspace=swanlab_workspace,
            reinit=True,
            name=trial_name,
            config=config,
        )
        swanlab_init_kwargs.update(self.kwargs)

        # 启动日志记录Actor
        self._start_logging_actor(trial, exclude_results, **swanlab_init_kwargs)

    def _start_logging_actor(
        self, trial: "Trial", exclude_results: List[str], **swanlab_init_kwargs
    ):
        """
        启动日志记录Actor
        
        为指定试验创建并启动后台日志记录Actor
        
        调用逻辑：
        log_trial_start -> _start_logging_actor -> 创建远程Actor -> Actor.run()
        
        Args:
            trial: 试验对象
            exclude_results: 要排除的结果键列表
            **swanlab_init_kwargs: SwanLab初始化参数
        """
        # 如果Actor已存在，重用（试验重启时可能发生）
        if trial in self._trial_logging_futures:
            return

        # 创建远程Actor类（如果还没有）
        if not self._remote_logger_class:
            env_vars = {}
            self._remote_logger_class = ray.remote(
                num_cpus=0,
                **_force_on_current_node(),
                runtime_env={"env_vars": env_vars},
                max_restarts=-1,
                max_task_retries=-1,
            )(self._logger_actor_cls)

        # 创建队列用于与Actor通信
        self._trial_queues[trial] = Queue(
            actor_options={
                "num_cpus": 0,
                **_force_on_current_node(),
                "max_restarts": -1,
                "max_task_retries": -1,
            }
        )
        
        # 创建Actor实例
        self._trial_logging_actors[trial] = self._remote_logger_class.remote(
            logdir=trial.local_path,
            queue=self._trial_queues[trial],
            exclude=exclude_results,
            to_config=self.AUTO_CONFIG_KEYS,
            **swanlab_init_kwargs,
        )
        
        # 启动Actor运行
        logging_future = self._trial_logging_actors[trial].run.remote()
        self._trial_logging_futures[trial] = logging_future
        self._logging_future_to_trial[logging_future] = trial

    def _signal_logging_actor_stop(self, trial: "Trial"):
        """
        向日志记录Actor发送停止信号
        
        Args:
            trial: 试验对象
        """
        self._trial_queues[trial].put((_QueueItem.END, None))

    def log_trial_result(self, iteration: int, trial: "Trial", result: Dict):
        """
        记录试验结果
        
        当试验产生新结果时调用，负责：
        1. 确保Actor已启动
        2. 清理结果数据
        3. 发送数据到Actor队列
        
        调用逻辑：
        Ray Tune -> log_trial_result -> 发送数据到队列 -> Actor._handle_result
        
        Args:
            iteration: 迭代次数
            trial: 试验对象
            result: 结果字典
        """
        # 如果Actor还没有启动，先启动
        if trial not in self._trial_logging_actors:
            self.log_trial_start(trial)

        # 清理结果数据
        result = _clean_log(result)
        # 发送到Actor队列
        self._trial_queues[trial].put((_QueueItem.RESULT, result))

    def log_trial_save(self, trial: "Trial"):
        """
        记录试验保存事件
        
        当试验保存检查点时调用
        
        Args:
            trial: 试验对象
        """
        if self.upload_checkpoints and trial.checkpoint:
            checkpoint_root = None
            if isinstance(trial.checkpoint.filesystem, pyarrow.fs.LocalFileSystem):
                checkpoint_root = trial.checkpoint.path

            if checkpoint_root:
                self._trial_queues[trial].put((_QueueItem.CHECKPOINT, checkpoint_root))

    def log_trial_end(self, trial: "Trial", failed: bool = False):
        """
        记录试验结束
        
        当试验结束时调用，负责：
        1. 向Actor发送停止信号
        2. 清理Actor资源
        
        调用逻辑：
        Ray Tune -> log_trial_end -> _signal_logging_actor_stop -> _cleanup_logging_actors
        
        Args:
            trial: 试验对象
            failed: 是否失败
        """
        self._signal_logging_actor_stop(trial=trial)
        self._cleanup_logging_actors()

    def _cleanup_logging_actor(self, trial: "Trial"):
        """
        清理单个试验的日志记录Actor
        
        Args:
            trial: 试验对象
        """
        del self._trial_queues[trial]
        del self._trial_logging_futures[trial]
        ray.kill(self._trial_logging_actors[trial])
        del self._trial_logging_actors[trial]

    def _cleanup_logging_actors(self, timeout: int = 0, kill_on_timeout: bool = False):
        """
        清理所有日志记录Actor
        
        等待Actor完成并清理资源
        
        Args:
            timeout: 等待超时时间
            kill_on_timeout: 超时时是否强制杀死Actor
        """
        futures = list(self._trial_logging_futures.values())
        done, remaining = ray.wait(futures, num_returns=len(futures), timeout=timeout)
        
        # 清理已完成的Actor
        for ready_future in done:
            finished_trial = self._logging_future_to_trial.pop(ready_future)
            self._cleanup_logging_actor(finished_trial)

        # 如果超时且需要强制杀死，清理剩余的Actor
        if kill_on_timeout:
            for remaining_future in remaining:
                trial = self._logging_future_to_trial.pop(remaining_future)
                self._cleanup_logging_actor(trial)

    def on_experiment_end(self, trials: List["Trial"], **info):
        """
        实验结束时的回调
        
        等待所有Actor完成对swanlab.finish的调用。
        这包括上传所有日志和工件到SwanLab。
        
        调用逻辑：
        Ray Tune -> on_experiment_end -> _cleanup_logging_actors
        
        Args:
            trials: 所有试验列表
            **info: 额外信息
        """
        self._cleanup_logging_actors(timeout=self._upload_timeout, kill_on_timeout=True)

    def __del__(self):
        """
        析构函数
        
        确保在对象销毁时清理所有资源
        """
        if ray.is_initialized():
            # 向所有Actor发送停止信号
            for trial in list(self._trial_logging_actors):
                self._signal_logging_actor_stop(trial=trial)

            # 清理所有Actor
            self._cleanup_logging_actors(timeout=2, kill_on_timeout=True)

        # 清空所有字典
        self._trial_logging_actors = {}
        self._trial_logging_futures = {}
        self._logging_future_to_trial = {}
        self._trial_queues = {}
