"""
@author: mx
@file: __init__.py
@time: 2026/10/28
@description: 直方图处理模块
"""

import math
from numbers import Integral
from typing import List, Optional, Tuple

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramValue
from swanlab.sdk.internal.context import TransformHistogram
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.transforms.histogram import (
    HistogramBinCount,
    HistogramBins,
    HistogramInput,
    NumericSequence,
)

MAX_BUCKETS = 512
MIN_BUCKETS = 1
MAX_COUNT = (1 << 64) - 1


class Histogram(TransformHistogram):
    """直方图，对一批样本的分段计数摘要。

    两种输入方式互斥：原始序列 sequence 由 SDK 内部分桶，预分桶结果 np_histogram 直接采用，
    提供 np_histogram 时 sequence、bins、num_bins 均被忽略。

    :param sequence: 原始样本序列，支持 list/tuple 与 numpy 数组。PyTorch Tensor 请先调用
        tensor.detach().cpu().numpy() 转换；带梯度或位于 GPU 上的 Tensor 不能直接传入。
    :param bins: 桶数（1..512 的 Python / NumPy 整数，按序列最小/最大值等宽分桶）
        或桶边界序列（N+1 个单调递增值，支持 numpy 数组），默认 64 个等宽桶
    :param num_bins: 向下兼容 wandb 的桶数参数，语义与整数形式的 bins 相同，与 bins 同时
        显式给出视为矛盾输入
    :param np_histogram: 预分桶结果 (counts, bin_edges)，提供时优先于 sequence
    """

    def __init__(
        self,
        sequence: Optional[HistogramInput] = None,
        bins: Optional[HistogramBins] = None,
        *,
        num_bins: Optional[HistogramBinCount] = None,
        np_histogram: Optional[Tuple[NumericSequence, NumericSequence]] = None,
    ):
        super().__init__()
        # 套娃加载：内层 Histogram 实例提供已计算的结果；内层结果是最终形态，无法重新分桶，
        # 因此外层再给任何分桶参数都视为矛盾输入
        inner = self._unwrap(sequence)
        if inner:
            if np_histogram is not None or bins is not None or num_bins is not None:
                raise ValueError("Cannot specify bins or np_histogram when nesting a Histogram instance")
            self.values: List[int] = inner["values"]
            self.bins: List[float] = inner["bins"]
            return

        if np_histogram is not None:
            self._init_from_np_histogram(np_histogram)
            return

        if sequence is None:
            raise ValueError("Either sequence or np_histogram must be provided")

        if num_bins is not None and (isinstance(num_bins, bool) or not isinstance(num_bins, Integral)):
            raise TypeError(f"num_bins must be an integer bucket count, got {type(num_bins).__name__}")
        if bins is not None and num_bins is not None:
            raise ValueError("bins and num_bins are mutually exclusive, got both")
        effective_bins = bins if bins is not None else (num_bins if num_bins is not None else 64)
        self._init_from_sequence(sequence, effective_bins)

    # ----------------------------------
    # 两种输入的计算与校验
    # ----------------------------------

    def _init_from_np_histogram(self, np_histogram) -> None:
        """预分桶输入：校验后原样采用，不调用 numpy 分桶。"""
        if not (isinstance(np_histogram, (tuple, list)) and len(np_histogram) == 2):
            raise ValueError("Expected np_histogram to be a tuple of (values, bin_edges)")
        values = self._validate_counts(np_histogram[0])
        bins = self._validate_edges(np_histogram[1])
        if len(bins) != len(values) + 1:
            raise ValueError(f"len(bins) must be len(values) + 1, got {len(bins)} edges and {len(values)} buckets")
        self.values = values
        self.bins = bins

    def _init_from_sequence(self, sequence, bins) -> None:
        """原始序列输入：过滤非有限值后调用 numpy 分桶。"""
        if isinstance(bins, bool):
            raise TypeError("bins must be an integer bucket count or a sequence of floats, got bool")

        if isinstance(bins, Integral):
            bins = int(bins)
            if not (MIN_BUCKETS <= bins <= MAX_BUCKETS):
                raise ValueError(f"The number of buckets must be between {MIN_BUCKETS} and {MAX_BUCKETS}, got {bins}")
            finite = self._finite_values(sequence)
            if len(finite) == 0:
                raise ValueError("sequence contains no finite values")
            counts, edges = vendor.np.histogram(finite, bins=bins)
            self.values = counts.tolist()
            self.bins = edges.tolist()
            return

        # 显式边界：先完成边界校验，再分桶并提示范围外丢弃的样本数
        edges_list = self._validate_edges(bins)
        finite = self._finite_values(sequence)
        counts, edges = vendor.np.histogram(finite, bins=edges_list)
        in_range = ((finite >= edges_list[0]) & (finite <= edges_list[-1])).sum()
        if in_range < len(finite):
            console.warning(
                f"{len(finite) - int(in_range)} value(s) in sequence fall outside the bins and are excluded"
            )
        self.values = counts.tolist()
        self.bins = edges.tolist()

    @staticmethod
    def _numeric_array(data, name: str) -> "vendor.np.ndarray":
        """校验一维实数输入，保留原始数值类型。"""
        if isinstance(data, (str, bytes)):
            raise TypeError(f"{name} must contain numeric values")
        arr = vendor.np.asarray(data)
        if arr.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional, got ndim={arr.ndim}")
        if arr.dtype.kind not in "biuf":
            raise TypeError(f"{name} must contain real numeric values")
        return arr

    @classmethod
    def _finite_values(cls, sequence) -> "vendor.np.ndarray":
        """丢弃原始序列中的 NaN/Inf 并提示丢弃数量。"""
        arr = cls._numeric_array(sequence, "sequence")
        finite = arr[vendor.np.isfinite(arr)]
        dropped = len(arr) - len(finite)
        if dropped > 0:
            console.warning(f"{dropped} non-finite value(s) in sequence are excluded from the histogram")
        return finite

    @staticmethod
    def _validate_counts(counts) -> List[int]:
        """校验每桶计数：一维实数、有限、非负整数，且不超过 UInt64 上限。"""
        if isinstance(counts, (str, bytes)):
            raise TypeError("Bucket counts must be numeric")
        # object 数组避免 NumPy 将混合整数输入提升为 float64，保留 UInt64 精度。
        arr = vendor.np.asarray(counts, dtype=object)
        if arr.ndim != 1:
            raise ValueError(f"Bucket counts must be one-dimensional, got ndim={arr.ndim}")
        values: List[int] = []
        for value in arr:
            if not isinstance(value, (int, float, vendor.np.integer, vendor.np.floating)):
                raise TypeError("Bucket counts must contain real numeric values")
            if not isinstance(value, (int, vendor.np.integer)) and not vendor.np.isfinite(value):
                raise ValueError("Bucket counts must not contain NaN or Inf")
            if value < 0:
                raise ValueError("Bucket counts must be non-negative")
            count = int(value)
            if value != count:
                raise ValueError("Bucket counts must be integers")
            if count > MAX_COUNT:
                raise ValueError(f"Bucket counts must not exceed UInt64 maximum {MAX_COUNT}")
            values.append(count)
        return values

    @classmethod
    def _validate_edges(cls, edges) -> List[float]:
        """校验桶边界：一维、数值、有限、单调递增，桶数范围 1..512，重复值警告。"""
        arr = cls._numeric_array(edges, "bins").astype(float)
        if not vendor.np.isfinite(arr).all():
            raise ValueError("bins must not contain NaN or Inf")
        edges_list = arr.tolist()
        if len(edges_list) < MIN_BUCKETS + 1:
            raise ValueError(f"At least {MIN_BUCKETS} bucket is required, got {len(edges_list) - 1} buckets")
        if len(edges_list) - 1 > MAX_BUCKETS:
            raise ValueError(f"The number of buckets must not exceed {MAX_BUCKETS}")
        diffs = [b - a for a, b in zip(edges_list, edges_list[1:])]
        if any(d < 0 for d in diffs):
            raise ValueError("bins must be monotonically increasing")
        if any(math.isclose(d, 0.0) for d in diffs):
            console.warning("bins contains duplicate values, zero-width buckets will be produced")
        return edges_list

    # ----------------------------------
    # TransformData 接口
    # ----------------------------------

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_HISTOGRAM

    def transform(self) -> HistogramValue:
        """将已计算的分桶结果转换为 Protobuf 格式。"""
        return HistogramValue(values=self.values, bins=self.bins)
