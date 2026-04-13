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
    get_swanlab_version,
    has_run,
    init,
    log,
    log_audio,
    log_image,
    log_text,
    log_video,
    login,
    logout,
    merge_settings,
    verify,
)

from . import utils

__version__ = get_swanlab_version()

__all__ = [
    # cmd
    "merge_settings",
    "init",
    "finish",
    "login",
    "logout",
    "verify",
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
