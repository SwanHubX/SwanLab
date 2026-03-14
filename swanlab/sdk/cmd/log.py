"""
@author: cunyue
@file: log.py
@time: 2026/3/14
@description: SwanLab SDK logging methods
"""

from typing import Any, Mapping, Optional, Union

from swanlab.sdk.cmd.helper import with_cmd_lock, with_run
from swanlab.sdk.internal.run import get_run
from swanlab.sdk.internal.run.data.transforms import Text


@with_cmd_lock
@with_run("log")
def log(data: Mapping[str, Any], step: Optional[int] = None):
    """Log metrics and data to the current run.

    :param data: Dictionary of metric names and values to log.

    :param step: Optional step number. If not provided, auto-increments.

    :raises RuntimeError: If called without an active run.

    Examples:

        Log single metric:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Log multiple metrics:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5, "accuracy": 0.95})
        >>> swanlab.finish()

        Log with explicit step:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5}, step=10)
        >>> swanlab.finish()
    """
    run = get_run()
    run.log(data, step)


@with_cmd_lock
@with_run("log_text")
def log_text(key: str, data: Union[str, Text], caption: Optional[str] = None, step: Optional[int] = None):
    """Log text data to the current run.

    :param key: The key for the text data.

    :param data: The text data itself or a Text object.

    :param caption: Optional caption for the text data.

    :param step: Optional step number. If not provided, auto-increments.

    :raises RuntimeError: If called without an active run.

    Examples:

        Log simple text:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log_text("output", "Training started")
        >>> swanlab.finish()

        Log text with caption:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log_text("prediction", "cat", caption="Model output")
        >>> swanlab.finish()

        Log text with step:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log_text("status", "Epoch 5 complete", step=5)
        >>> swanlab.finish()
    """
    run = get_run()
    run.log_text(key, data, caption, step)
