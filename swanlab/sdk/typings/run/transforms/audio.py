"""
@author: cunyue
@file: audio.py
@time: 2026/3/17 19:35
@description: 音频处理模块
"""

from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab import vendor
    from swanlab.sdk.internal.run.transforms.audio import Audio

AudioDataType = Union[
    "Audio",
    str,
    "vendor.np.ndarray",
]

AudioDatasType = Union[
    "Audio",
    List["Audio"],
    str,
    List[str],
    "vendor.np.ndarray",
    List["vendor.np.ndarray"],
]


AudioRateType = int

AudioRatesType = Union[int, List[int]]
