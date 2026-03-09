from swanlab.sdk import Settings, finish, init, login, merge_settings, utils
from swanlab.sdk.utils.version import get_swanlab_version

__all__ = [
    "merge_settings",
    "Settings",
    "init",
    "finish",
    "login",
    "utils",
]


__version__ = get_swanlab_version()
