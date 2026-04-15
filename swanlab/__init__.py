from swanlab.sdk import (
    Audio,
    Callback,
    Image,
    Run,
    Settings,
    Text,
    Video,
    async_log,
    config,
    define_scalar,
    finish,
    get_run,
    has_run,
    helper,
    init,
    log,
    log_audio,
    log_image,
    log_text,
    log_video,
    login,
    merge_settings,
)

from . import utils

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
    "Image",
    "Video",
    # utils
    "utils",
    "Settings",
    "Callback",
]


def __getattr__(name: str):
    if name == "run":
        try:
            return get_run()
        except RuntimeError:
            return None
    raise AttributeError(f"module 'swanlab' has no attribute {name!r}")
