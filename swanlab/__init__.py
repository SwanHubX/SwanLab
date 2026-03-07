from swanlab.sdk import Settings, finish, init, merge_settings
from swanlab.sdk.pkg.version import get_swanlab_version

__all__ = [
    "merge_settings",
    "Settings",
    "init",
    "finish",
]


__version__ = get_swanlab_version()
