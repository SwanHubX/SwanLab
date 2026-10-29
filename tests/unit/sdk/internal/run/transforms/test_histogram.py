"""
@author: mx
@file: test_histogram.py
@time: 2026/10/28
@description: 测试直方图转换模块 Histogram，覆盖需求文档边界情况表
"""

from unittest.mock import patch

import numpy as np
import pytest

from swanlab import Histogram
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord

# ---------------------------------- 原始序列输入 ----------------------------------


class TestSequenceInput:
    def test_default_bins(self):
        """默认 64 个等宽桶"""
        h = Histogram(np.arange(128, dtype=float))
        assert len(h.values) == 64
        assert len(h.bins) == 65
        assert sum(h.values) == 128

    def test_int_bins(self):
        """指定桶数"""
        h = Histogram([1, 2, 3, 4, 5], bins=4)
        assert len(h.values) == 4
        assert len(h.bins) == 5

    @pytest.mark.parametrize("bins", [[-1.5, 0.0, 0.75, 3.0], (-1.5, 0.0, 0.75, 3.0)])
    def test_explicit_edges(self, bins):
        """显式边界序列，含小数"""
        h = Histogram([0.5, 1.5, 2.5], bins=bins)
        assert h.values == [0, 1, 2]
        assert h.bins == [-1.5, 0.0, 0.75, 3.0]

    def test_numpy_semantics_last_bin_closed(self):
        """除最后一桶外左闭右开，最后一桶包含右端点（numpy 语义）"""
        h = Histogram([0.0, 1.0, 2.0, 3.0], bins=[0, 1, 2, 3])
        assert h.values == [1, 1, 2]

    def test_num_bins_compat(self):
        """num_bins 向下兼容 wandb 的桶数参数"""
        h = Histogram([1, 2, 3], num_bins=2)
        assert len(h.values) == 2

    def test_tensor_sequence(self):
        """torch tensor 序列（仅当 torch 已安装）"""
        torch = pytest.importorskip("torch")
        h = Histogram(torch.tensor([1.0, 2.0, 3.0]), bins=2)
        assert sum(h.values) == 3

    def test_bins_and_num_bins_conflict(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            Histogram([1, 2, 3], bins=4, num_bins=4)

    def test_bucket_count_out_of_range(self):
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2, 3], bins=513)
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2, 3], bins=0)

    def test_edges_bucket_count_out_of_range(self):
        edges = [float(i) for i in range(514)]
        with pytest.raises(ValueError, match="must not exceed 512"):
            Histogram([1, 2, 3], bins=edges)

    def test_empty_edges(self):
        with pytest.raises(ValueError, match="At least 1 bucket"):
            Histogram([1, 2, 3], bins=[1.0])

    @pytest.mark.parametrize("bins", ["64", True, np.array([0.0, 1.0, 2.0])])
    def test_bins_wrong_type(self, bins):
        with pytest.raises(TypeError):
            Histogram([1, 2, 3], bins=bins)

    @pytest.mark.parametrize("num_bins", ["64", 2.0, True])
    def test_num_bins_wrong_type(self, num_bins):
        with pytest.raises(TypeError, match="num_bins must be an int"):
            Histogram([1, 2, 3], num_bins=num_bins)

    def test_sequence_missing(self):
        with pytest.raises(ValueError, match="Either sequence or np_histogram"):
            Histogram()

    @pytest.mark.parametrize(
        "sequence",
        ["abc", ["a", "b", "c"], ["1", "2"], [1, "2"], np.array([1 + 2j])],
    )
    def test_non_numeric_sequence(self, sequence):
        with pytest.raises(TypeError):
            Histogram(sequence)

    @pytest.mark.parametrize("sequence", [[[1, 2], [3, 4]]])
    def test_two_dimensional_sequence(self, sequence):
        with pytest.raises(ValueError, match="one-dimensional"):
            Histogram(sequence)


class TestSequenceNonFinite:
    """非有限值属于非人为情况：丢弃并提示，不报错"""

    def test_nan_dropped(self):
        h = Histogram(np.array([1.0, np.nan, 3.0]), bins=2)
        assert sum(h.values) == 2

    def test_inf_dropped(self):
        h = Histogram(np.array([1.0, np.inf, 3.0]), bins=2)
        assert sum(h.values) == 2

    def test_all_non_finite(self):
        with pytest.raises(ValueError, match="no finite values"):
            Histogram(np.array([np.nan, np.inf]))

    def test_all_non_finite_with_explicit_edges(self):
        with patch("swanlab.sdk.internal.run.transforms.histogram.console.warning") as warning:
            h = Histogram([np.nan, np.inf, -np.inf], bins=[0.0, 1.0, 2.0])

        assert h.values == [0, 0]
        assert h.bins == [0.0, 1.0, 2.0]
        assert list(h.transform().values) == [0, 0]
        warning.assert_called_once_with("3 non-finite value(s) in sequence are excluded from the histogram")

    def test_empty_sequence_with_explicit_edges(self):
        h = Histogram([], bins=[0.0, 1.0, 2.0])
        assert h.values == [0, 0]

    def test_values_outside_edges_dropped(self):
        """显式边界之外的样本静默丢弃（与 numpy 一致）"""
        h = Histogram([0.5, 1.5, 2.5, 9.0], bins=[0.0, 1.0, 2.0])
        assert h.values == [1, 1]

    def test_all_values_outside_edges(self):
        """全部在边界外：全 0 直方图，允许"""
        h = Histogram([10.0, 20.0], bins=[0.0, 1.0, 2.0])
        assert h.values == [0, 0]


# ---------------------------------- 预分桶输入 ----------------------------------


class TestNpHistogramInput:
    def test_precomputed(self):
        counts, edges = np.histogram([1, 2, 3, 4], bins=2)
        h = Histogram(np_histogram=(counts, edges))
        assert h.values == counts.tolist()
        assert h.bins == edges.tolist()

    @pytest.mark.parametrize("np_histogram", [([1, 2], [0, 1], [2, 3]), (), ([1],)])
    def test_not_a_pair(self, np_histogram):
        with pytest.raises(ValueError, match="Expected np_histogram"):
            Histogram(np_histogram=np_histogram)

    def test_length_mismatch(self):
        with pytest.raises(ValueError, match="len\\(bins\\) must be len\\(values\\) \\+ 1"):
            Histogram(np_histogram=([1, 2], [0, 1]))

    def test_non_monotonic_edges(self):
        with pytest.raises(ValueError, match="monotonically increasing"):
            Histogram(np_histogram=([1, 1], [0.0, 2.0, 1.0]))

    def test_nan_edges(self):
        with pytest.raises(ValueError, match="must not contain NaN or Inf"):
            Histogram(np_histogram=([1, 1], [0.0, np.nan, 2.0]))

    def test_inf_edges(self):
        with pytest.raises(ValueError, match="must not contain NaN or Inf"):
            Histogram(np_histogram=([1, 1], [0.0, np.inf, 2.0]))

    def test_negative_counts(self):
        with pytest.raises(ValueError, match="non-negative"):
            Histogram(np_histogram=([-1, 1], [0.0, 1.0, 2.0]))

    def test_non_integer_counts(self):
        with pytest.raises(ValueError, match="must be integers"):
            Histogram(np_histogram=([0.5, 1], [0.0, 1.0, 2.0]))

    def test_integer_valued_float_counts(self):
        """整数值的 float 计数接受并转为 int"""
        h = Histogram(np_histogram=([1.0, 2.0], [0.0, 1.0, 2.0]))
        assert h.values == [1, 2]

    @pytest.mark.parametrize("counts", [["a"], ["2"], [1, "2"], np.array([1 + 2j])])
    def test_non_numeric_counts(self, counts):
        with pytest.raises(TypeError):
            Histogram(np_histogram=(counts, [0.0, 1.0]))

    @pytest.mark.parametrize("edges", [["0", "1"], [0, "1"], np.array([0 + 1j, 1 + 2j])])
    def test_non_numeric_edges(self, edges):
        with pytest.raises(TypeError):
            Histogram(np_histogram=([1], edges))

    @pytest.mark.parametrize("counts", [[[1, 2]]])
    def test_two_dimensional_counts(self, counts):
        with pytest.raises(ValueError, match="one-dimensional"):
            Histogram(np_histogram=(counts, [0.0, 1.0, 2.0]))

    @pytest.mark.parametrize(
        "counts",
        [
            [2**53 + 1],
            [0, 2**64 - 1],
            np.array([2**64 - 1], dtype=np.uint64),
            [1.0, np.uint64(2**64 - 1)],
        ],
    )
    def test_uint64_counts_preserve_precision(self, counts):
        expected = [int(count) for count in counts]
        h = Histogram(np_histogram=(counts, list(range(len(counts) + 1))))

        assert h.values == expected
        value = h.transform()
        restored = type(value).FromString(value.SerializeToString())
        assert list(restored.values) == expected

    @pytest.mark.parametrize("count", [2**64, 2**64 + 1, float(2**64), 10**400])
    def test_counts_outside_uint64_range(self, count):
        with pytest.raises(ValueError, match="UInt64 maximum"):
            Histogram(np_histogram=([count], [0.0, 1.0]))

    @pytest.mark.parametrize("count", [np.nan, np.inf, -np.inf])
    def test_non_finite_counts(self, count):
        with pytest.raises(ValueError, match="must not contain NaN or Inf"):
            Histogram(np_histogram=([count], [0.0, 1.0]))

    def test_empty_counts(self):
        with pytest.raises(ValueError, match="At least 1 bucket"):
            Histogram(np_histogram=([], [1.0]))

    def test_duplicate_edges_allowed(self):
        """重复边界产生零宽桶：允许"""
        h = Histogram(np_histogram=([1, 0], [0.0, 1.0, 1.0]))
        assert h.bins == [0.0, 1.0, 1.0]


# ---------------------------------- 优先级与套娃 ----------------------------------


class TestPriority:
    def test_np_histogram_overrides_sequence(self):
        h = Histogram([1, 2, 3], np_histogram=([5], [0, 1]))
        assert h.values == [5]

    def test_np_histogram_ignores_bins(self):
        h = Histogram([1, 2, 3], bins=99, np_histogram=([5], [0, 1]))
        assert h.values == [5]

    def test_nesting_without_args(self):
        inner = Histogram([1, 2, 3], bins=2)
        h = Histogram(inner)
        assert h.values == inner.values
        assert h.bins == inner.bins

    def test_nesting_with_bins_conflict(self):
        inner = Histogram([1, 2, 3], bins=2)
        with pytest.raises(ValueError, match="nesting"):
            Histogram(inner, bins=3)


# ---------------------------------- Protobuf 转换 ----------------------------------


class TestTransform:
    def test_column_type(self):
        assert Histogram.column_type() == ColumnType.COLUMN_TYPE_HISTOGRAM

    def test_transform(self):
        h = Histogram([1, 2, 3], bins=2)
        value = h.transform()
        assert list(value.values) == h.values
        assert list(value.bins) == h.bins

    def test_build_data_record(self):
        from google.protobuf.timestamp_pb2 import Timestamp

        ts = Timestamp()
        ts.GetCurrentTime()
        h = Histogram([1, 2, 3], bins=2)
        record = Histogram.build_data_record(key="train/w", step=3, timestamp=ts, data=h.transform())
        assert isinstance(record, HistogramRecord)
        assert record.key == "train/w"
        assert record.step == 3
        assert record.type == ColumnType.COLUMN_TYPE_HISTOGRAM
        assert list(record.value.values) == h.values
