"""
@author: cunyue
@file: __init__.py
@time: 2026/3/10 13:37
@description: SwanLab 第三方库集成
部分第三方库并非SwanLab必须依赖，在引入时需要判断是否已安装
考虑到可维护性和性能，我们采用延迟导入的方式，仅在实际使用时才导入第三方库
"""
# pyright: reportMissingImports=false

import importlib
from typing import TYPE_CHECKING, Any

# 1. Type hinting block: Only executed by static type checkers (e.g., Pyright, MyPy, IDEs)
if TYPE_CHECKING:
    import accelerate
    import accelerate.tracking
    import boto3
    import fastai
    import fastai.learner
    import fastcore
    import fastcore.basics
    import imageio
    import keras
    import lightning
    import lightning.pytorch
    import lightning.pytorch.loggers
    import lightning.pytorch.utilities
    import matplotlib
    import matplotlib.figure
    import mlflow
    import mmengine
    import mmengine.config
    import mmengine.registry
    import mmengine.visualization.vis_backend
    import moviepy
    import numpy as np
    import paddlenlp
    import paddlenlp.trainer.trainer
    import pandas as pd
    import PIL
    import PIL.Image
    import pynvml
    import rdkit
    import rdkit.Chem
    import rdkit.Chem.AllChem
    import sklearn
    import sklearn.metrics
    import soundfile
    import stable_baselines3
    import stable_baselines3.common
    import swanboard
    import tensorboard
    import tensorboard.backend.event_processing
    import tensorboard.util
    import tensorboardX
    import torch
    import torchtune
    import torchtune.utils.metric_logging
    import torchvision
    import transformers
    import ultralytics
    import wandb
    import xgboost
    import xgboost.callback


# 2. Expose the available modules for IDE auto-completion
__all__ = [
    "imageio",
    "matplotlib",
    "moviepy",
    "np",
    "PIL",
    "rdkit",
    "sklearn",
    "soundfile",
    "swanboard",
    "tensorboard",
    "tensorboardX",
    "boto3",
    "torch",
    "torchvision",
    "pynvml",
    "wandb",
    "mlflow",
    # these are extra dependencies which are not in [project.optional-dependencies]
    "pd",
    # framework integrations
    "accelerate",
    "fastai",
    "fastcore",
    "keras",
    "lightning",
    "mmengine",
    "paddlenlp",
    "transformers",
    "ultralytics",
    "xgboost",
    "stable_baselines3",
    "torchtune",
]

# 3. Lazy import mapping: Actual module paths
_LAZY_IMPORTS = {
    "imageio": "imageio",
    "matplotlib": "matplotlib",
    "moviepy": "moviepy",
    "np": "numpy",
    "PIL": "PIL",
    "rdkit": "rdkit",
    "sklearn": "sklearn",
    "soundfile": "soundfile",
    "swanboard": "swanboard",
    "tensorboard": "tensorboard",
    "tensorboardX": "tensorboardX",
    "boto3": "boto3",
    "torch": "torch",
    "torchvision": "torchvision",
    "pynvml": "pynvml",
    # these are extra dependencies which are not in [project.optional-dependencies]
    "pd": "pandas",
    # framework integrations — users install these themselves
    "accelerate": "accelerate",
    "fastai": "fastai",
    "fastcore": "fastcore",
    "keras": "keras",
    "lightning": "lightning",
    "mmengine": "mmengine",
    "paddlenlp": "paddlenlp",
    "transformers": "transformers",
    "ultralytics": "ultralytics",
    "stable_baselines3": "stable_baselines3",
    "torchtune": "torchtune",
    "wandb": "wandb",
    "mlflow": "mlflow",
    "xgboost": "xgboost",
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

# 5. Submodule imports: some packages require submodules to be imported explicitly
# so their attributes are accessible (e.g. PIL.Image must be imported for PIL.Image to work)
_SUBMODULE_IMPORTS = {
    "PIL": ["PIL.Image"],
    "matplotlib": ["matplotlib.figure"],
    "sklearn": ["sklearn.metrics"],
    "rdkit": ["rdkit.Chem", "rdkit.Chem.AllChem"],
    "accelerate": ["accelerate.tracking"],
    "fastai": ["fastai.learner", "fastai.callback.hook"],
    "fastcore": ["fastcore.basics"],
    "keras": ["keras.callbacks"],
    "lightning": ["lightning.pytorch", "lightning.pytorch.loggers", "lightning.pytorch.utilities"],
    "mmengine": ["mmengine.config", "mmengine.registry", "mmengine.visualization.vis_backend"],
    "transformers": ["transformers.trainer_callback"],
    "xgboost": ["xgboost.callback"],
    "paddlenlp": ["paddlenlp.trainer.trainer"],
    "stable_baselines3": ["stable_baselines3.common"],
    "torchtune": ["torchtune.utils.metric_logging"],
    "tensorboard": [
        "tensorboard.backend.event_processing",
        "tensorboard.backend.event_processing.event_file_loader",
        "tensorboard.util",
    ],
}


# 6. Module-level __getattr__ for lazy loading (PEP 562)
def _install_hint(name: str, module_path: str) -> str:
    """顶层包缺失时的安装提示文案。"""
    extra_tag = _EXTRA_DEPS.get(name)
    if extra_tag:
        return (
            f"The '{name}' feature requires additional dependencies. "
            f"To enable it, please install the '{extra_tag}' extra by running:\n"
            f'    pip install "swanlab[{extra_tag}]"'
        )
    return (
        f"The '{name}' feature requires the '{module_path}' package, "
        f"which is not currently installed. Please install it by running:\n"
        f"    pip install {module_path}"
    )


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
        except ImportError as e:
            # 仅当顶层包本身缺失时才提示安装，否则报告真实错误
            if isinstance(e, ModuleNotFoundError) and e.name == module_path:
                raise ImportError(_install_hint(name, module_path)) from e
            raise ImportError(f"The '{name}' feature failed to import '{module_path}': {e}") from e

        # Import required submodules so their attributes are accessible on the parent package
        for submodule_path in _SUBMODULE_IMPORTS.get(name, []):
            try:
                importlib.import_module(submodule_path)
            except ImportError as e:
                raise ImportError(f"The '{name}' feature requires '{submodule_path}' of '{module_path}': {e}") from e

        # Cache the imported object in the module's global namespace
        globals()[name] = obj
        return obj

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
