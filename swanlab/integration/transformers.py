"""
Docs: https://docs.swanlab.cn/zh/guide_cloud/integration/integration-huggingface-transformers.html
"""

from typing import Optional, List, Dict, Union, Any
import swanlab

try:
    from transformers.trainer_callback import TrainerCallback
except ImportError:
    raise RuntimeError(
        "This contrib module requires Transformers to be installed. "
        "Please install it with command: \n pip install transformers"
    )


def rewrite_logs(d):
    new_d = {}
    eval_prefix = "eval_"
    eval_prefix_len = len(eval_prefix)
    test_prefix = "test_"
    test_prefix_len = len(test_prefix)
    for k, v in d.items():
        if k.startswith(eval_prefix):
            new_d["eval/" + k[eval_prefix_len:]] = v
        elif k.startswith(test_prefix):
            new_d["test/" + k[test_prefix_len:]] = v
        else:
            new_d["train/" + k] = v
    return new_d


class SwanLabCallback(TrainerCallback):
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[str] = None,
        **kwargs: Any,
    ):
        self._initialized = False
        self._experiment = swanlab

        self._swanlab_init: Dict[str, Any] = {
            "project": project,
            "workspace": workspace,
            "experiment_name": experiment_name,
            "description": description,
            "logdir": logdir,
            "mode": mode,
        }

        self._swanlab_init.update(**kwargs)

        self._project = self._swanlab_init.get("project")
        self._workspace = self._swanlab_init.get("workspace")
        self._experiment_name = self._swanlab_init.get("experiment_name")
        self._description = self._swanlab_init.get("decsription")
        self._logdir = self._swanlab_init.get("logdir")
        self._mode = self._swanlab_init.get("mode")

    def setup(self, args, state, model, **kwargs):
        self._initialized = True

        if not state.is_world_process_zero:
            return

        # 如果没有注册过实验
        if self._experiment.get_run() is None:
            self._experiment.init(**self._swanlab_init)

        if args:
            combined_dict = {**args.to_sanitized_dict()}
            for key, value in combined_dict.items():
                self._experiment.config.set(key, value)

    def on_train_begin(self, args, state, control, model=None, **kwargs):
        if not self._initialized:
            self.setup(args, state, model, **kwargs)

    def on_train_end(self, args, state, control, model=None, tokenizer=None, **kwargs):
        pass

    def on_log(self, args, state, control, model=None, logs=None, **kwargs):
        if not self._initialized:
            self.setup(args, state, model, **kwargs)

        if state.is_world_process_zero:
            logs = rewrite_logs(logs)
            self._experiment.log(logs, step=state.global_step)
