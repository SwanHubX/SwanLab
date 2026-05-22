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
    merge_callbacks,
    merge_settings,
    pkg,
    plot,
    save,
    sync,
)

from . import utils
from .api import Api
from .deprecated.callbacks import register_callbacks
from .deprecated.echarts import confusion_matrix, pr_curve, roc_curve

__version__ = pkg.helper.get_swanlab_version()

__all__ = [
    # cmd
    "merge_callbacks",
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
    "save",
    "async_log",
    "sync",
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
    "register_callbacks",
]


def __getattr__(name: str):
    if name == "run":
        try:
            return get_run()
        except RuntimeError:
            return None
    raise AttributeError(f"module 'swanlab' has no attribute {name!r}")
