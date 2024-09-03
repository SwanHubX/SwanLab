"""
Docs: https://docs.swanlab.cn/zh/guide_cloud/integration/integration-huggingface-transformers.html
"""
import warnings
from .transformers import *

# 只显示一次DeprecationWarning
warnings.simplefilter('once', DeprecationWarning)

# 发出弃用警告
warnings.warn(
    "The module 'huggingface' is deprecated and will be removed in future versions. "
    "Please update your imports to use 'transformers' instead.",
    DeprecationWarning,
    stacklevel=2,
)
