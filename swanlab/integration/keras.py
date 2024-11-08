import sys
from typing import Any, Dict, Optional, Union
import tensorflow as tf
try:
    from tensorflow.keras.callbacks import Callback
    from tensorflow.keras.utils import model_to_dot
except ImportError as exc:
    try:
        from keras.callbacks import Callback
        from keras.utils import model_to_dot
    except ImportError:
        msg = """
        keras package not found.

        As Keras is now part of Tensorflow you should install it by running
            pip install tensorflow"""
        raise ModuleNotFoundError(msg) from exc

import swanlab

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

LogStrategy = Literal["epoch", "batch"]


class SwanLabLogger(Callback):
    def __init__(
        self,
        log_freq: Union[LogStrategy, int] = "epoch",
        initial_global_step: int = 0,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)

        if swanlab.get_run() is None:
            raise RuntimeError(
                "You must call `swanlab.init()` before SwanLabLogger()"
            )

        if log_freq == "batch":
            log_freq = 1

        self.logging_batch_wise = isinstance(log_freq, int)
        self.log_freq: Any = log_freq if self.logging_batch_wise else None
        self.global_batch = 0
        self.global_step = initial_global_step

    def _get_lr(self) -> Union[float, None]:
        if isinstance(
            self.model.optimizer.learning_rate,
            (tf.Variable, tf.Tensor),
        ) or (
            hasattr(self.model.optimizer.learning_rate, "shape")
            and self.model.optimizer.learning_rate.shape == ()
        ):
            return float(self.model.optimizer.learning_rate.numpy().item())
        try:
            return float(
                self.model.optimizer.learning_rate(step=self.global_step).numpy().item()
            )
        except Exception as e:
            return None

    def on_epoch_end(self, epoch: int, logs: Optional[Dict[str, Any]] = None) -> None:
        """在每个epoch结束时调用"""
        logs = dict() if logs is None else {f"epoch/{k}": v for k, v in logs.items()}

        logs["epoch/epoch"] = epoch

        lr = self._get_lr()
        if lr is not None:
            logs["epoch/learning_rate"] = lr

        swanlab.log(logs)

    def on_batch_end(self, batch: int, logs: Optional[Dict[str, Any]] = None) -> None:
        self.global_step += 1
        """与 `on_train_batch_end` 的向后兼容别名"""
        if self.logging_batch_wise and batch % self.log_freq == 0:
            logs = {f"batch/{k}": v for k, v in logs.items()} if logs else {}
            logs["batch/batch_step"] = self.global_batch

            lr = self._get_lr()
            if lr is not None:
                logs["batch/learning_rate"] = lr

            swanlab.log(logs)

            self.global_batch += self.log_freq

    def on_train_batch_end(
        self, batch: int, logs: Optional[Dict[str, Any]] = None
    ) -> None:
        """在 `fit` 方法的训练batch结束时调用"""
        self.on_batch_end(batch, logs if logs else {})
