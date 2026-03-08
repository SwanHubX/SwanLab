from swanlab.sdk import Settings, finish, init, login, merge_settings
from swanlab.sdk.pkg.version import get_swanlab_version

__all__ = [
    "merge_settings",
    "Settings",
    "init",
    "finish",
    "login",
]


__version__ = get_swanlab_version()
