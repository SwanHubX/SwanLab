"""
@author: mx
@file: __init__.py
@time: 2026/10/28
@description: 直方图处理模块
"""

import math
from numbers import Integral
from typing import Any, List, Optional, Tuple

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramValue
from swanlab.sdk.internal.context import TransformHistogram
from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.typings.run.transforms.histogram import (
    HistogramBins,
    HistogramInput,
    NumericSequence,
)


def _is_torch_tensor(obj: Any) -> bool:
    """通过类型名检测 PyTorch Tensor/Parameter，避免强制导入 torch。"""
    typename = f"{obj.__class__.__module__}.{obj.__class__.__name__}"
    return typename.startswith("torch.") and any(name in typename for name in ("Tensor", "Variable", "Parameter"))


def _to_numpy(data: Any) -> Any:
    """将 PyTorch Tensor 转为 numpy 数组，兼容 GPU、requires_grad 与 bfloat16。"""
    if not _is_torch_tensor(data):
        return data
    tensor = data.detach().cpu()
    dtype_name = str(getattr(tensor, "dtype", ""))
    # numpy 不支持 bfloat16 / float8；直方图边界最终为 double，转 float32 不影响分桶语义
    if "bfloat16" in dtype_name or "float8" in dtype_name:
        tensor = tensor.float()
    return tensor.numpy()


class Histogram(TransformHistogram):
    """Histogram of binned sample counts.

    The two input modes are mutually exclusive: ``sequence`` is binned by the
    SDK, while a precomputed ``np_histogram`` is used as-is. When
    ``np_histogram`` is provided, ``sequence``, ``bins`` and ``num_bins`` are
    ignored.

    Parameters
    ----------
    sequence: Sequence[float], numpy.ndarray, or Histogram, optional
        Raw samples to bin. PyTorch tensors are converted automatically with
        ``detach().cpu()``, including GPU and requires_grad tensors.
    bins: int or Sequence[float], optional
        Number of buckets (1..512) for equal-width binning, or explicit bucket
        edges (N+1 monotonically increasing values). ``None`` means 64 buckets.
    np_histogram: tuple[Sequence[int], Sequence[float]], optional
        Precomputed ``(counts, bin_edges)``; takes precedence over ``sequence``.
    **kwargs
        ``num_bins`` is accepted as a wandb-compatible alias for an integer
        ``bins``. It is used only when ``bins`` is omitted; an explicitly
        provided ``bins`` wins and emits a warning when the two differ. Other
        unknown keyword arguments are ignored.

    Notes
    -----
    With explicit ``bins``, an empty sequence or a sequence containing only
    non-finite values yields an all-zero histogram. Integer ``bins`` needs a
    finite range, so it raises ``ValueError`` when no finite values remain.

    Type and shape errors (unsupported objects, non-numeric data, wrong
    dimensionality) are logged and mark the instance invalid instead of
    raising; logging an invalid histogram skips it. Value errors such as an
    out-of-range bucket count or non-monotonic edges still raise. Missing
    optional dependencies (e.g. ``numpy``) raise ``ImportError`` with an
    install hint instead of being treated as invalid input.
    """

    DEFAULT_BINS = 64
    MIN_BUCKETS = 1
    MAX_BUCKETS = 512
    MAX_COUNT = (1 << 64) - 1

    def __init__(
        self,
        sequence: Optional[HistogramInput] = None,
        bins: Optional[HistogramBins] = None,
        *,
        np_histogram: Optional[Tuple[NumericSequence, NumericSequence]] = None,
        **kwargs: Any,
    ):
        super().__init__()
        self.values: List[int] = []
        self.bins: List[float] = []
        self._error: Optional[BaseException] = None

        # wandb 兼容：num_bins 是 bins 的别名，其余未知 kwargs 静默忽略
        num_bins = kwargs.pop("num_bins", None)

        # 参数组合错误属于调用方式问题，不受类型兜底影响，直接抛出
        if isinstance(sequence, Histogram):
            if np_histogram is not None or bins is not None or num_bins is not None:
                raise ValueError("Cannot specify bins or np_histogram when nesting a Histogram instance")
        elif np_histogram is None and sequence is None:
            raise ValueError("Either sequence or np_histogram must be provided")

        # 无效类型/形状 log 时跳过；只兜输入校验错误，
        # 可选依赖缺失（如 numpy）必须上抛，否则安装提示会被当成无效输入吞掉
        with safe.block(
            TypeError,
            ValueError,
            message="Invalid histogram input, the metric will be ignored",
            on_error=self._capture_error,
            write_to_tty=False,
        ):
            prepared = self._prepare_inputs(sequence, bins, num_bins, np_histogram)

        if self._error is not None:
            return

        self._build(prepared)

    def _capture_error(self, error: BaseException) -> None:
        """记录类型/形状校验失败原因，供 log 阶段跳过该实例。"""
        self._error = error

    @property
    def is_valid(self) -> bool:
        """输入类型与形状是否通过校验。"""
        return self._error is None

    # ----------------------------------
    # 输入归一化与数值校验
    # ----------------------------------

    @classmethod
    def _prepare_inputs(cls, sequence, bins, num_bins, np_histogram):
        """只做类型与形状归一化，不做数值语义校验。

        本阶段抛出的异常由构造函数中的 safe.block 捕获并转为 invalid 状态。
        """
        if isinstance(sequence, Histogram):
            if not sequence.is_valid:
                raise ValueError("Nested Histogram instance is invalid")
            return ("nested", sequence)

        if np_histogram is not None:
            if not (isinstance(np_histogram, (tuple, list)) and len(np_histogram) == 2):
                raise ValueError("Expected np_histogram to be a tuple of (values, bin_edges)")
            return ("np", cls._prepare_counts(np_histogram[0]), cls._prepare_edges(np_histogram[1]))

        # wandb 兼容：bins 未指定时 num_bins 生效；bins 已指定时以 bins 为准
        if num_bins is not None:
            if bins is None:
                if isinstance(num_bins, bool) or not isinstance(num_bins, Integral):
                    raise TypeError(f"num_bins must be an integer bucket count, got {type(num_bins).__name__}")
                bins = int(num_bins)
            else:
                console.warning(f"Conflicting bins={bins!r} and num_bins={num_bins!r}, using bins={bins!r}")

        if isinstance(bins, bool):
            raise TypeError("bins must be an integer bucket count or a sequence of floats, got bool")
        if isinstance(bins, Integral):
            prepared_bins = int(bins)
        elif bins is None:
            prepared_bins = cls.DEFAULT_BINS
        else:
            prepared_bins = cls._prepare_edges(bins)

        return ("seq", cls._numeric_array(sequence, "sequence"), prepared_bins)

    @staticmethod
    def _prepare_counts(counts) -> "vendor.np.ndarray":
        """校验 counts 的类型与维度，保留原始数值类型。"""
        if isinstance(counts, (str, bytes)):
            raise TypeError("Bucket counts must be numeric")
        # object 数组避免 NumPy 将混合整数输入提升为 float64，保留 UInt64 精度。
        arr = vendor.np.asarray(_to_numpy(counts), dtype=object)
        if arr.ndim != 1:
            raise ValueError(f"Bucket counts must be one-dimensional, got ndim={arr.ndim}")
        for value in arr:
            if not isinstance(value, (int, float, vendor.np.integer, vendor.np.floating)):
                raise TypeError("Bucket counts must contain real numeric values")
        return arr

    @classmethod
    def _prepare_edges(cls, edges) -> "vendor.np.ndarray":
        """校验 bins 的类型与维度，并统一为 float 数组。"""
        return cls._numeric_array(edges, "bins").astype(float)

    def _build(self, prepared) -> None:
        """执行数值语义校验与分桶计算；这里的 ValueError 仍然抛给调用方。"""
        mode = prepared[0]
        if mode == "nested":
            inner: Histogram = prepared[1]
            self.values = list(inner.values)
            self.bins = list(inner.bins)
            return

        if mode == "np":
            values = self._validate_counts(prepared[1])
            bins = self._validate_edges(prepared[2])
            if len(bins) != len(values) + 1:
                raise ValueError(f"len(bins) must be len(values) + 1, got {len(bins)} edges and {len(values)} buckets")
            self.values = values
            self.bins = bins
            return

        sequence: "vendor.np.ndarray" = prepared[1]
        bins = prepared[2]
        if isinstance(bins, int):
            if not (self.MIN_BUCKETS <= bins <= self.MAX_BUCKETS):
                raise ValueError(
                    f"The number of buckets must be between {self.MIN_BUCKETS} and {self.MAX_BUCKETS}, got {bins}"
                )
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
        arr = vendor.np.asarray(_to_numpy(data))
        if arr.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional, got ndim={arr.ndim}")
        if arr.dtype.kind not in "biuf":
            raise TypeError(f"{name} must contain real numeric values")
        return arr

    @classmethod
    def _finite_values(cls, arr: "vendor.np.ndarray") -> "vendor.np.ndarray":
        """丢弃原始序列中的 NaN/Inf 并提示丢弃数量。"""
        finite = arr[vendor.np.isfinite(arr)]
        dropped = len(arr) - len(finite)
        if dropped > 0:
            console.warning(f"{dropped} non-finite value(s) in sequence are excluded from the histogram")
        return finite

    @classmethod
    def _validate_counts(cls, arr: "vendor.np.ndarray") -> List[int]:
        """校验已归一化的每桶计数：有限、非负整数，且不超过 UInt64 上限。"""
        values: List[int] = []
        for value in arr:
            if not isinstance(value, (int, vendor.np.integer)) and not vendor.np.isfinite(value):
                raise ValueError("Bucket counts must not contain NaN or Inf")
            if value < 0:
                raise ValueError("Bucket counts must be non-negative")
            count = int(value)
            if value != count:
                raise ValueError("Bucket counts must be integers")
            if count > cls.MAX_COUNT:
                raise ValueError(f"Bucket counts must not exceed UInt64 maximum {cls.MAX_COUNT}")
            values.append(count)
        return values

    @classmethod
    def _validate_edges(cls, arr: "vendor.np.ndarray") -> List[float]:
        """校验已归一化的桶边界：有限、单调递增，桶数范围 1..512，重复值警告。"""
        if not vendor.np.isfinite(arr).all():
            raise ValueError("bins must not contain NaN or Inf")
        edges_list = arr.tolist()
        if len(edges_list) < cls.MIN_BUCKETS + 1:
            raise ValueError(f"At least {cls.MIN_BUCKETS} bucket is required, got {len(edges_list) - 1} buckets")
        if len(edges_list) - 1 > cls.MAX_BUCKETS:
            raise ValueError(f"The number of buckets must not exceed {cls.MAX_BUCKETS}")
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
        if self._error is not None:
            raise ValueError(f"Invalid histogram input: {self._error}") from self._error
        return HistogramValue(values=self.values, bins=self.bins)
