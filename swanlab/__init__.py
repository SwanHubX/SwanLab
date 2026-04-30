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
    # utils
    "utils",
    "echarts",
    "Settings",
    "Callback",
    # Api
    "Api",
]


def __getattr__(name: str):
    if name == "run":
        try:
            return get_run()
        except RuntimeError:
            return None
    raise AttributeError(f"module 'swanlab' has no attribute {name!r}")
