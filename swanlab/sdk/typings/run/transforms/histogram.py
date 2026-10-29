"""
@author: mx
@file: histogram.py
@time: 2026/10/28
@description: 直方图接受参数类型定义
"""

from typing import TYPE_CHECKING, Sequence, Union

if TYPE_CHECKING:
    from swanlab import vendor
    from swanlab.sdk.internal.run.transforms.histogram import Histogram

Number = Union[int, float]
NumericSequence = Union[Sequence[Number], "vendor.np.ndarray"]
HistogramInput = Union[NumericSequence, "Histogram"]
HistogramBinCount = Union[int, "vendor.np.integer"]
HistogramBins = Union[HistogramBinCount, NumericSequence]
