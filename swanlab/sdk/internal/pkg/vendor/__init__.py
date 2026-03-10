"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 13:37
@description: SwanLab 第三方库集成
部分第三方库并非SwanLab必须依赖，在引入时需要判断是否已安装
考虑到可维护性和性能，我们采用延迟导入的方式，仅在实际使用时才导入第三方库
"""

import importlib
from typing import TYPE_CHECKING, Any

# 1. Type hinting block: Only executed by static type checkers (e.g., Pyright, MyPy, IDEs)
if TYPE_CHECKING:
    import boto3
    import imageio
    import matplotlib
    import moviepy
    import numpy as np
    import pandas as pd
    import PIL
    import rdkit
    import soundfile
    import swanboard


# 2. Expose the available modules for IDE auto-completion
__all__ = [
    "imageio",
    "matplotlib",
    "moviepy",
    "np",
    "PIL",
    "rdkit",
    "soundfile",
    "swanboard",
    "boto3",
    # these are extra dependencies which are not in [project.optional-dependencies]
    "pd",
]

# 3. Lazy import mapping: Actual module paths
_LAZY_IMPORTS = {
    "imageio": "imageio",
    "matplotlib": "matplotlib",
    "moviepy": "moviepy",
    "np": "numpy",
    "PIL": "PIL",
    "rdkit": "rdkit",
    "soundfile": "soundfile",
    "swanboard": "swanboard",
    "boto3": "boto3",
    # these are extra dependencies which are not in [project.optional-dependencies]
    "pd": "pandas",
}

# 4. Optional dependencies mapping: Maps imported names to SwanLab's 'extras'
# This is strictly based on the [project.optional-dependencies] in pyproject.toml
_EXTRA_DEPS = {
    # [project.optional-dependencies.media]
    "soundfile": "media",
    "PIL": "media",  # Package 'pillow' is imported as 'PIL'
    "matplotlib": "media",
    "np": "media",
    "moviepy": "media",
    "imageio": "media",
    "rdkit": "media",
    # [project.optional-dependencies.dashboard]
    "swanboard": "dashboard",
    # [project.optional-dependencies.s3]
    "boto3": "s3",
}


# 5. Module-level __getattr__ for lazy loading (PEP 562)
def __getattr__(name: str) -> Any:
    if name in _LAZY_IMPORTS:
        module_path = _LAZY_IMPORTS[name]

        try:
            # Handle relative imports for internal integration modules
            if module_path.startswith("."):
                module = importlib.import_module(module_path, package=__name__)
                obj = getattr(module, name)
            else:
                # Handle direct third-party library imports
                obj = importlib.import_module(module_path)

            # Cache the imported object in the module's global namespace
            globals()[name] = obj
            return obj

        except ImportError as e:
            extra_tag = _EXTRA_DEPS.get(name)

            if extra_tag:
                error_msg = (
                    f"The '{name}' feature requires additional dependencies. "
                    f"To enable it, please install the '{extra_tag}' extra by running:\n"
                    f'    pip install "swanlab[{extra_tag}]"'
                )
            else:
                # Fallback: if not mapped in _EXTRA_DEPS, suggest the underlying package
                underlying_pkg = module_path.strip(".")
                error_msg = (
                    f"The '{name}' feature requires the '{underlying_pkg}' package, "
                    f"which is not currently installed. Please install it by running:\n"
                    f"    pip install {underlying_pkg}"
                )

            raise ImportError(error_msg) from e

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
