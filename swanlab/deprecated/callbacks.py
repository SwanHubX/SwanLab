"""
@author: cunyue
@description: callbacks相关的即将被弃用的命令
"""

import warnings

from typing_extensions import deprecated

from swanlab.sdk.cmd.merge_callbacks import merge_callbacks


@deprecated("use `swanlab.merge_callbacks()` instead")
def register_callbacks(*args, **kwargs):
    warnings.warn(
        "`swanlab.register_callbacks()` is deprecated, use `swanlab.merge_callbacks()` instead",
        FutureWarning,
        stacklevel=2,
    )
    return merge_callbacks(*args, **kwargs)
