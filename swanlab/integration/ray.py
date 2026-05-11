"""
Docs: https://docs.swanlab.cn/guide_cloud/integration/integration-ray.html

Usage:
------
from ray import tune
from swanlab.integration.ray import SwanLabLoggerCallback

cb = SwanLabLoggerCallback(project="Ray_Project")
tuner = tune.Tuner(
    train_func,
    param_space={"lr": tune.grid_search([0.001, 0.01, 0.1, 1.0]), "epochs": 10},
    run_config=tune.RunConfig(callbacks=[cb]),
)
results = tuner.fit()
------
"""

from __future__ import annotations

import enum
import os
import warnings
from numbers import Number
from typing import Any, Dict, List, Optional, Tuple
from urllib.error import HTTPError

import swanlab
import swanlab.vendor
from swanlab import Callback

try:
    # Ray must be imported explicitly for submodules in 2.x (lazy loading)
    import ray  # noqa: F401
    import ray.air  # noqa: F401
    import ray.train  # noqa: F401
    import ray.tune  # noqa: F401
    import ray.tune.logger  # noqa: F401
    import ray.util  # noqa: F401
    from ray.air.constants import TRAINING_ITERATION as _TRAINING_ITERATION
    from ray.air.util.node import _force_on_current_node
    from ray.train._internal.syncer import DEFAULT_SYNC_TIMEOUT as _DEFAULT_SYNC_TIMEOUT
    from ray.tune.experiment import Trial as _Trial
    from ray.tune.logger import LoggerCallback as _LoggerCallback
    from ray.tune.utils import flatten_dict as _flatten_dict
    from ray.util.queue import Queue

    _ray_logger = ray.tune.logger
except ImportError as e:
    raise ImportError(
        "The Ray integration requires the 'ray' package. Please install it by running:\n    pip install 'ray[tune]'"
    ) from e

np = swanlab.vendor.np

# Environment constants
SWANLAB_API_KEY = "SWANLAB_API_KEY"
SWANLAB_MODE = "SWANLAB_MODE"
SWANLAB_PROJ_NAME = "SWANLAB_PROJ_NAME"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_swanlab_project(project: Optional[str] = None) -> Optional[str]:
    if not project and os.environ.get(SWANLAB_PROJ_NAME):
        project = os.environ.get(SWANLAB_PROJ_NAME)
    return project


def _is_allowed_type(obj: Any) -> bool:
    if isinstance(obj, np.ndarray) and obj.size == 1:
        return isinstance(obj.item(), Number)
    return isinstance(obj, Number)


def _clean_log(obj: Any) -> Any:
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
        return int(fallback)
    except ValueError:
        pass
    return fallback


# ---------------------------------------------------------------------------
# Queue items
# ---------------------------------------------------------------------------


class _QueueItem(enum.Enum):
    END = enum.auto()
    RESULT = enum.auto()
    CHECKPOINT = enum.auto()


# ---------------------------------------------------------------------------
# Ray Actor
# ---------------------------------------------------------------------------


class _SwanLabLoggingActor:
    def __init__(
        self,
        logdir: str,
        queue: Queue,
        exclude: List[str],
        to_config: List[str],
        *args: Any,
        **kwargs: Any,
    ) -> None:
        self._swanlab = swanlab
        os.chdir(logdir)
        self.queue = queue
        self._exclude = set(exclude)
        self._to_config = set(to_config)
        self.args = args
        self.kwargs = kwargs
        self._trial_name = self.kwargs.get("name", "unknown")
        self._logdir = logdir

    def run(self) -> None:
        run = self._swanlab.init(*self.args, **self.kwargs)
        run.config.trial_log_path = self._logdir

        while True:
            item_type, item_content = self.queue.get()
            if item_type == _QueueItem.END:
                break
            if item_type == _QueueItem.CHECKPOINT:
                continue

            log, config_update = self._handle_result(item_content)
            try:
                self._swanlab.config.update(config_update, allow_val_change=True)
                self._swanlab.log(log, step=log.get(_TRAINING_ITERATION))
            except HTTPError as e:
                _ray_logger.warning("Failed to log result to swanlab: {}".format(str(e)))  # type: ignore

        self._swanlab.finish()

    def _handle_result(self, result: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        config_update = result.get("config", {}).copy()
        log: Dict[str, Any] = {}
        flat_result = _flatten_dict(result, delimiter="/")
        for k, v in flat_result.items():
            if any(k.startswith(item + "/") or k == item for item in self._exclude):
                continue
            elif any(k.startswith(item + "/") or k == item for item in self._to_config):
                config_update[k] = v
            elif not _is_allowed_type(v):
                continue
            else:
                log[k] = v
        config_update.pop("callbacks", None)
        return log, config_update


# ---------------------------------------------------------------------------
# SwanLabLoggerCallback
# ---------------------------------------------------------------------------


class SwanLabLoggerCallback(Callback, _LoggerCallback):
    _exclude_results = ["done", "should_checkpoint"]
    AUTO_CONFIG_KEYS = [
        "trial_id",
        "experiment_tag",
        "node_ip",
        "experiment_id",
        "hostname",
        "pid",
        "date",
    ]
    _logger_actor_cls = _SwanLabLoggingActor

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
        upload_timeout: int = _DEFAULT_SYNC_TIMEOUT,
        **kwargs: Any,
    ) -> None:
        if save_checkpoints:
            warnings.warn(
                "`save_checkpoints` is deprecated. Use `upload_checkpoints` instead.",
                DeprecationWarning,
            )
            upload_checkpoints = save_checkpoints

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

        self._remote_logger_class: Any = None
        self._trial_logging_actors: Dict[_Trial, Any] = {}
        self._trial_logging_futures: Dict[_Trial, Any] = {}
        self._logging_future_to_trial: Dict[Any, _Trial] = {}
        self._trial_queues: Dict[_Trial, Queue] = {}

    # --- swanlab.Callback ---

    @property
    def name(self) -> str:
        return "swanlab-integration-ray"

    def on_run_initialized(self, run_dir: Any, path: str, **kwargs: Any) -> None:
        run = _get_active_run()
        if run is not None:
            run.config["FRAMEWORK"] = "ray"

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs: Any) -> None:
        pass

    # --- Ray LoggerCallback ---

    def setup(self, *args: Any, **kwargs: Any) -> None:
        self.project = _get_swanlab_project(self.project)
        if not self.project:
            raise ValueError(
                f"Please pass the project name as argument or through the {SWANLAB_PROJ_NAME} environment variable."
            )

    def log_trial_start(self, trial: _Trial) -> None:
        config = trial.config.copy()
        config.pop("callbacks", None)
        exclude_results = self._exclude_results.copy()
        exclude_results += self.excludes
        if not self.log_config:
            exclude_results += ["config"]

        trial_id = trial.trial_id if trial else None
        trial_name = str(trial) if trial else None

        swanlab_init_kwargs = dict(
            project=self.project,
            workspace=self.workspace,
            name=trial_name,
            config=_clean_log(config),
        )
        swanlab_init_kwargs.update(self.kwargs)
        self._start_logging_actor(trial, exclude_results, **swanlab_init_kwargs)

    def _start_logging_actor(self, trial: _Trial, exclude_results: List[str], **swanlab_init_kwargs: Any) -> None:
        if trial in self._trial_logging_futures:
            return
        if not self._remote_logger_class:
            env_vars = {}
            self._remote_logger_class = ray.remote(
                num_cpus=0,
                **_force_on_current_node(),
                runtime_env={"env_vars": env_vars},
                max_restarts=-1,
                max_task_retries=-1,
            )(self._logger_actor_cls)

        self._trial_queues[trial] = Queue(
            actor_options={
                "num_cpus": 0,
                **_force_on_current_node(),
                "max_restarts": -1,
                "max_task_retries": -1,
            }
        )
        self._trial_logging_actors[trial] = self._remote_logger_class.remote(
            logdir=trial.local_path,
            queue=self._trial_queues[trial],
            exclude=exclude_results,
            to_config=self.AUTO_CONFIG_KEYS,
            **swanlab_init_kwargs,
        )
        logging_future = self._trial_logging_actors[trial].run.remote()
        self._trial_logging_futures[trial] = logging_future
        self._logging_future_to_trial[logging_future] = trial

    def _signal_logging_actor_stop(self, trial: _Trial) -> None:
        self._trial_queues[trial].put((_QueueItem.END, None))

    def log_trial_result(self, iteration: int, trial: _Trial, result: Dict[str, Any]) -> None:
        if trial not in self._trial_logging_actors:
            self.log_trial_start(trial)
        result = _clean_log(result)
        self._trial_queues[trial].put((_QueueItem.RESULT, result))

    def log_trial_save(self, trial: _Trial) -> None:
        if self.upload_checkpoints and trial.checkpoint:
            try:
                import pyarrow.fs

                if isinstance(trial.checkpoint.filesystem, pyarrow.fs.LocalFileSystem):
                    checkpoint_root = trial.checkpoint.path
                    if checkpoint_root:
                        self._trial_queues[trial].put((_QueueItem.CHECKPOINT, checkpoint_root))
            except ImportError:
                pass

    def log_trial_end(self, trial: _Trial, failed: bool = False) -> None:
        self._signal_logging_actor_stop(trial)
        self._cleanup_logging_actors()

    def _cleanup_logging_actor(self, trial: _Trial) -> None:
        del self._trial_queues[trial]
        del self._trial_logging_futures[trial]
        ray.kill(self._trial_logging_actors[trial])
        del self._trial_logging_actors[trial]

    def _cleanup_logging_actors(self, timeout: int = 0, kill_on_timeout: bool = False) -> None:
        futures = list(self._trial_logging_futures.values())
        done, remaining = ray.wait(futures, num_returns=len(futures), timeout=timeout)
        for ready_future in done:
            finished_trial = self._logging_future_to_trial.pop(ready_future)
            self._cleanup_logging_actor(finished_trial)
        if kill_on_timeout:
            for remaining_future in remaining:
                trial = self._logging_future_to_trial.pop(remaining_future)
                self._cleanup_logging_actor(trial)

    def on_experiment_end(self, trials: List[_Trial], **info: Any) -> None:
        self._cleanup_logging_actors(timeout=self._upload_timeout, kill_on_timeout=True)

    def __del__(self) -> None:
        try:
            if ray.is_initialized():
                for trial in list(self._trial_logging_actors):
                    self._signal_logging_actor_stop(trial)
                self._cleanup_logging_actors(timeout=2, kill_on_timeout=True)
            self._trial_logging_actors = {}
            self._trial_logging_futures = {}
            self._logging_future_to_trial = {}
            self._trial_queues = {}
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------


def _get_active_run():
    try:
        return swanlab.get_run()
    except RuntimeError:
        return None


__all__ = ["SwanLabLoggerCallback"]
