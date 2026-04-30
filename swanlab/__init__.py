import warnings

from typing_extensions import deprecated

from swanlab.sdk import (
    Audio,
    Callback,
    ECharts,
    Image,
    Molecule,
    Object3D,
    Run,
    Settings,
    Text,
    Video,
    async_log,
    config,
    define_scalar,
    echarts,
    finish,
    get_run,
    has_run,
    helper,
    init,
    log,
    log_audio,
    log_echarts,
    log_image,
    log_molecule,
    log_object3d,
    log_text,
    log_video,
    login,
    merge_settings,
    plot,
)

from . import utils
from .api import Api

__version__ = helper.get_swanlab_version()

__all__ = [
    # cmd
    "merge_settings",
    "init",
    "finish",
    "login",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "log_echarts",
    "log_object3d",
    "log_molecule",
    "define_scalar",
    "async_log",
    # run
    "run",  # type: ignore [no-redef]
    "Run",
    "has_run",
    "get_run",
    # config
    "config",
    # data
    "Text",
    "Audio",
    "ECharts",
    "Image",
    "Molecule",
    "Object3D",
    "Video",
    "plot",
    "echarts",
    # utils
    "utils",
    "Settings",
    "Callback",
    # Api
    "Api",
    # deprecated
    "roc_curve",
    "pr_curve",
    "confusion_matrix",
]


def __getattr__(name: str):
    if name == "run":
        try:
            return get_run()
        except RuntimeError:
            return None
    raise AttributeError(f"module 'swanlab' has no attribute {name!r}")


# ---------------------------------- deprecated ----------------------------------


@deprecated("use `swanlab.echarts.roc_curve()` instead")
def roc_curve(*args, **kwargs):
    warnings.warn(
        "`swanlab.roc_curve()` is deprecated, use `swanlab.echarts.roc_curve()` instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return echarts.roc_curve(*args, **kwargs)


@deprecated("use `swanlab.echarts.pr_curve()` instead")
def pr_curve(*args, **kwargs):
    warnings.warn(
        "`swanlab.pr_curve()` is deprecated, use `swanlab.echarts.pr_curve()` instead", DeprecationWarning, stacklevel=2
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
