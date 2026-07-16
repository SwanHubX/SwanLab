"""
@author: cunyue
@description: 部分即将弃用的echarts语法糖
"""

import warnings

from typing_extensions import deprecated

from swanlab.sdk.internal.run.transforms import echarts


@deprecated("use `swanlab.echarts.roc_curve()` instead")
def roc_curve(*args, **kwargs):
    warnings.warn(
        "`swanlab.roc_curve()` is deprecated, use `swanlab.echarts.roc_curve()` instead",
        FutureWarning,
        stacklevel=2,
    )
    return echarts.roc_curve(*args, **kwargs)


@deprecated("use `swanlab.echarts.pr_curve()` instead")
def pr_curve(*args, **kwargs):
    warnings.warn(
        "`swanlab.pr_curve()` is deprecated, use `swanlab.echarts.pr_curve()` instead", FutureWarning, stacklevel=2
    )
    return echarts.pr_curve(*args, **kwargs)


@deprecated("use `swanlab.echarts.confusion_matrix()` instead")
def confusion_matrix(*args, **kwargs):
    warnings.warn(
        "`swanlab.confusion_matrix()` is deprecated, use `swanlab.echarts.confusion_matrix()` instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return echarts.confusion_matrix(*args, **kwargs)
