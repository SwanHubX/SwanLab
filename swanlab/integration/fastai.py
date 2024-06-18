"""
Docs: https://docs.swanlab.cn/zh/guide_cloud/integration/integration-fastai.html
"""

try:
    from fastai.learner import Callback
    from fastcore.basics import store_attr, detuplify, ignore_exceptions
    from fastai.callback.hook import total_params
except ImportError:
    raise RuntimeError(
        "This module requires `fastai` to be installed. " "Please install it with command: \n pip install fastai"
    )

from typing import Optional, Any
import swanlab
from swanlab.log import swanlog as swl


class SwanLabCallback(Callback):
    def __init__(
        self,
        project: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        workspace: Optional[str] = None,
        config: Optional[dict] = None,
        mode: Optional[str] = None,
        logdir: Optional[str] = None,
        **kwargs: Any,
    ):
        store_attr()
        self._experiment = swanlab

        self.project = project
        self.experiment_name = experiment_name
        self.workspace = workspace
        self.config = config
        self.description = description
        self.mode = mode
        self.logdir = logdir
        self.train_suffix = "train"
        self.summary_suffix = "summary"

    def setup_swanlab(self):
        if self._experiment.get_run() is None:
            self._experiment.init(
                project=self.project,
                workspace=self.workspace,
                experiment_name=self.experiment_name,
                description=self.description,
                config=self.config,
                mode=self.mode,
                logdir=self.logdir,
            )

    def before_fit(self):
        # print("=" * 10 + "before_fit" + "=" * 10)
        self.setup_swanlab()
        configs_log = self.gather_args()
        formatted_config = SwanLabCallback.format_config(configs_log)
        self._experiment.config.update(formatted_config)
        self._swanlab_step = 0

    def after_batch(self):
        if self.training:
            self._swanlab_step += 1
            swanlab.log(
                {f"{self.train_suffix}/loss": self.loss.item(), f"{self.train_suffix}/step": self._swanlab_step},
                step=self._swanlab_step,
            )
            for i, h in enumerate(self.opt.hypers):
                for k, v in h.items():
                    swanlab.log({f"{self.train_suffix}/{k}_{i}": v}, self._swanlab_step)

    def before_epoch(self):
        # print("=" * 10 + "before_epoch" + "=" * 10)
        for metric in self.metrics:
            metric.reset()

    def after_epoch(self):
        # print("=" * 10 + "after_epoch" + "=" * 10)
        for name, value in zip(self.recorder.metric_names, self.recorder.log):
            if value is not None:
                swanlab.log({f"{self.summary_suffix}/{name}": value})

    def __del__(self):
        # 如果实验已经结束，且实验状态为0，即RUNNING状态，则关闭实验
        if self._experiment.Run.get_state().value == 0:
            swanlab.finish()

    def gather_args(self):
        "Gather config parameters accessible to the learner"
        cb_args = {f"{cb}": getattr(cb, "__stored_args__", True) for cb in self.cbs if cb != self}
        args = {"Learner": self.learn, **cb_args}
        try:
            n_inp = self.dls.train.n_inp
            args["n_inp"] = n_inp
            xb = self.dls.valid.one_batch()[:n_inp]
            args.update(
                {f"input {n+1} dim {i+1}": d for n in range(n_inp) for i, d in enumerate(list(detuplify(xb[n]).shape))}
            )
        except Exception:
            swl.warning("Failed to gather input dimensions")
        with ignore_exceptions():
            args["batch_size"] = self.dls.bs
            args["batch_per_epoch"] = len(self.dls.train)
            args["model_parameters"] = total_params(self.model)[0]
            args["device"] = self.dls.device.type
            args["frozen"] = bool(self.opt.frozen_idx)
            args["frozen_idx"] = self.opt.frozen_idx
            args["dataset/tfms"] = f"{self.dls.dataset.tfms}"
            args["dls/after_item"] = f"{self.dls.after_item}"
            args["dls/before_batch"] = f"{self.dls.before_batch}"
            args["dls/after_batch"] = f"{self.dls.after_batch}"
        return args

    @classmethod
    def format_config(cls, config):
        "Format config parameters for logging"
        for key, value in config.items():
            if isinstance(value, dict):
                config[key] = SwanLabCallback.format_config(value)
            else:
                config[key] = SwanLabCallback.format_config_value(value)
        return config

    @classmethod
    def format_config_value(cls, value):
        if isinstance(value, list):
            return [SwanLabCallback.format_config_value(item) for item in value]
        elif hasattr(value, "__stored_args__"):
            return {**SwanLabCallback.format_config(value.__stored_args__), "_name": value}
        return value
