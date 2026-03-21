"""
swanlab.core_python — 顶层快捷导入包。
将 swanlab.sdk.internal.core_python 的公开接口重新导出，
使外部可以通过 swanlab.core_python.uploader 等短路径访问。
"""

import importlib
import sys

# 将 swanlab.sdk.internal.core_python 的子模块注册为 swanlab.core_python 的子模块
_REAL_PREFIX = "swanlab.sdk.internal.core_python"
_SHORT_PREFIX = "swanlab.core_python"

# 需要暴露的子模块
_SUBMODULES = [
    "uploader",
    "uploader.batch",
    "uploader.model",
    "uploader.upload",
    "uploader.thread",
    "uploader.thread.log_collector",
    "uploader.thread.start_thread",
    "uploader.thread.utils",
]

for _sub in _SUBMODULES:
    _real_name = f"{_REAL_PREFIX}.{_sub}"
    _short_name = f"{_SHORT_PREFIX}.{_sub}"
    if _short_name not in sys.modules:
        try:
            _mod = importlib.import_module(_real_name)
            sys.modules[_short_name] = _mod
        except ImportError:
            pass
