"""
@author: mx
@file: test_histogram.py
@time: 2026/10/28
@description: 测试直方图转换模块 Histogram，覆盖需求文档边界情况表
"""

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from swanlab import Histogram, Text
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import HistogramRecord

# ---------------------------------- 原始序列输入 ----------------------------------


def _make_invalid(*args, **kwargs):
    """构造一个类型/形状非法的 Histogram，并断言它被标记为无效"""
    with patch("swanlab.sdk.internal.run.transforms.histogram.console.trace") as trace:
        h = Histogram(*args, **kwargs)
    assert not h.is_valid
    assert h.values == [] and h.bins == []
    trace.assert_called_once()
    return h


class TestSequenceInput:
    def test_default_bins(self):
        """默认 64 个等宽桶"""
        h = Histogram(np.arange(128, dtype=float))
        assert len(h.values) == 64
        assert len(h.bins) == 65
        assert sum(h.values) == 128

    @pytest.mark.parametrize("bins", [4, np.int64(4), np.uint64(4)])
    def test_int_bins(self, bins):
        """指定桶数"""
        h = Histogram([1, 2, 3, 4, 5], bins=bins)
        assert len(h.values) == 4
        assert len(h.bins) == 5

    @pytest.mark.parametrize("bins", [[-1.5, 0.0, 0.75, 3.0], (-1.5, 0.0, 0.75, 3.0), np.array([-1.5, 0.0, 0.75, 3.0])])
    def test_explicit_edges(self, bins):
        """显式边界序列，含小数"""
        h = Histogram([0.5, 1.5, 2.5], bins=bins)
        assert h.values == [0, 1, 2]
        assert h.bins == [-1.5, 0.0, 0.75, 3.0]

    def test_range_edges(self):
        h = Histogram([0.5, 1.5, 2.5], bins=range(4))
        assert h.values == [1, 1, 1]
        assert h.bins == [0.0, 1.0, 2.0, 3.0]

    def test_numpy_semantics_last_bin_closed(self):
        """除最后一桶外左闭右开，最后一桶包含右端点（numpy 语义）"""
        h = Histogram([0.0, 1.0, 2.0, 3.0], bins=[0, 1, 2, 3])
        assert h.values == [1, 1, 2]

    @pytest.mark.parametrize("num_bins", [2, np.int64(2), np.uint64(2)])
    def test_num_bins_compat(self, num_bins):
        """num_bins 向下兼容 wandb 的桶数参数"""
        h = Histogram([1, 2, 3], num_bins=num_bins)
        assert len(h.values) == 2

    def test_tensor_sequence(self):
        """torch tensor 序列（仅当 torch 已安装）"""
        torch = pytest.importorskip("torch")
        h = Histogram(torch.tensor([1.0, 2.0, 3.0]), bins=2)
        assert sum(h.values) == 3

    def test_tensor_sequence_requires_grad(self):
        """带梯度的 tensor 自动 detach"""
        torch = pytest.importorskip("torch")
        h = Histogram(torch.tensor([1.0, 2.0, 3.0], requires_grad=True), bins=2)
        assert sum(h.values) == 3

    def test_tensor_sequence_bfloat16(self):
        """bfloat16 tensor 转为 float32 后分桶"""
        torch = pytest.importorskip("torch")
        h = Histogram(torch.tensor([1.0, 2.0, 3.0], dtype=torch.bfloat16, requires_grad=True), bins=2)
        assert sum(h.values) == 3

    def test_bins_and_num_bins_conflict(self):
        """同时传入时 bins 优先，不抛异常"""
        with patch("swanlab.sdk.internal.run.transforms.histogram.console.warning") as warning:
            h = Histogram([1, 2, 3], bins=2, num_bins=4)
        assert len(h.values) == 2
        warning.assert_called_once()
        assert "using bins=2" in warning.call_args[0][0]

    def test_num_bins_used_when_bins_not_set(self):
        """bins 保持默认值时，num_bins 生效（wandb 兼容）"""
        h = Histogram([1, 2, 3], num_bins=2)
        assert len(h.values) == 2

    def test_explicit_default_bins_still_wins(self):
        """显式传 bins=64 时，bins 优先于 num_bins"""
        with patch("swanlab.sdk.internal.run.transforms.histogram.console.warning") as warning:
            h = Histogram([1, 2, 3], bins=64, num_bins=2)
        assert len(h.values) == 64
        warning.assert_called_once()
        assert "using bins=64" in warning.call_args[0][0]

    def test_unknown_kwargs_ignored(self):
        """其余 wandb 参数静默忽略"""
        h = Histogram([1, 2, 3], num_bins=2, foo="bar")
        assert len(h.values) == 2

    def test_bucket_count_out_of_range(self):
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2, 3], bins=513)
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2, 3], bins=0)

    @pytest.mark.parametrize("count", [np.int64(0), np.int64(-1), np.uint64(513), np.uint64(2**64 - 1)])
    @pytest.mark.parametrize("parameter", ["bins", "num_bins"])
    def test_numpy_bucket_count_out_of_range(self, count, parameter):
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2, 3], **{parameter: count})

    def test_edges_bucket_count_out_of_range(self):
        edges = [float(i) for i in range(514)]
        with pytest.raises(ValueError, match="must not exceed 512"):
            Histogram([1, 2, 3], bins=edges)

    def test_empty_edges(self):
        with pytest.raises(ValueError, match="At least 1 bucket"):
            Histogram([1, 2, 3], bins=[1.0])

    @pytest.mark.parametrize("bins", ["64", True, [0, "1"], np.array([0 + 1j, 1 + 2j])])
    def test_bins_wrong_type(self, bins):
        _make_invalid(sequence=[1, 2, 3], bins=bins)

    @pytest.mark.parametrize("bins", [np.array([[0, 1], [2, 3]]), np.array(2)])
    def test_edges_wrong_dimension(self, bins):
        _make_invalid(sequence=[1, 2, 3], bins=bins)

    @pytest.mark.parametrize("num_bins", ["64", 2.0, True])
    def test_num_bins_wrong_type(self, num_bins):
        _make_invalid(sequence=[1, 2, 3], num_bins=num_bins)

    def test_sequence_missing(self):
        with pytest.raises(ValueError, match="Either sequence or np_histogram"):
            Histogram()

    @pytest.mark.parametrize(
        "sequence",
        ["abc", ["a", "b", "c"], ["1", "2"], [1, "2"], np.array([1 + 2j])],
    )
    def test_non_numeric_sequence(self, sequence):
        _make_invalid(sequence)

    @pytest.mark.parametrize("sequence", [[[1, 2], [3, 4]]])
    def test_two_dimensional_sequence(self, sequence):
        _make_invalid(sequence)


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


# ---------------------------------- 类型校验兜底 ----------------------------------


class TestTypeValidationFailOpen:
    """类型/形状错误只记录日志并标记无效，不抛给训练线程"""

    @pytest.mark.parametrize(
        "kwargs",
        [
            dict(sequence="abc"),
            dict(sequence=["a", "b"]),
            dict(sequence=[[1, 2], [3, 4]]),
            dict(sequence=5),
            dict(sequence=[1, 2], bins="64"),
            dict(sequence=[1, 2], bins=True),
            dict(sequence=[1, 2], num_bins=2.0),
            dict(np_histogram=123),
            dict(np_histogram=(["1"], [0.0, 1.0])),
            dict(np_histogram=([1], ["0", "1"])),
        ],
    )
    def test_invalid_type_marks_invalid(self, kwargs):
        with patch("swanlab.sdk.internal.run.transforms.histogram.console.trace") as trace:
            h = Histogram(**kwargs)
        assert not h.is_valid
        assert h.values == []
        assert h.bins == []
        trace.assert_called_once()

    def test_semantic_errors_still_raise(self):
        with pytest.raises(ValueError, match="between 1 and 512"):
            Histogram([1, 2], bins=0)
        with pytest.raises(ValueError, match="monotonically increasing"):
            Histogram([1, 2], bins=[0.0, 2.0, 1.0])


class TestOptionalDependencyMissing:
    """可选依赖缺失不能被当作无效输入吞掉，安装提示必须抛给调用方"""

    def test_missing_numpy_propagates(self, monkeypatch):
        from swanlab.sdk.internal.run.transforms import histogram as histogram_module

        class _MissingNumpy:
            def __getattr__(self, name):
                raise ImportError("The 'np' feature requires additional dependencies")

        monkeypatch.setattr(histogram_module, "vendor", SimpleNamespace(np=_MissingNumpy()))

        with pytest.raises(ImportError, match="requires additional dependencies"):
            Histogram([1.0, 2.0])


# ---------------------------------- 预分桶输入 ----------------------------------


class TestNpHistogramInput:
    def test_precomputed(self):
        counts, edges = np.histogram([1, 2, 3, 4], bins=2)
        h = Histogram(np_histogram=(counts, edges))
        assert h.values == counts.tolist()
        assert h.bins == edges.tolist()

    def test_tensor_counts_and_edges(self):
        """np_histogram 的 counts / edges 也支持 torch tensor"""
        torch = pytest.importorskip("torch")
        counts = torch.tensor([2.0, 1.0], requires_grad=True)
        edges = torch.tensor([0.0, 1.0, 2.0])
        h = Histogram(np_histogram=(counts, edges))
        assert h.values == [2, 1]
        assert h.bins == [0.0, 1.0, 2.0]

    @pytest.mark.parametrize("np_histogram", [([1, 2], [0, 1], [2, 3]), (), ([1],)])
    def test_not_a_pair(self, np_histogram):
        _make_invalid(np_histogram=np_histogram)

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
        _make_invalid(np_histogram=(counts, [0.0, 1.0]))

    @pytest.mark.parametrize("edges", [["0", "1"], [0, "1"], np.array([0 + 1j, 1 + 2j])])
    def test_non_numeric_edges(self, edges):
        _make_invalid(np_histogram=([1], edges))

    @pytest.mark.parametrize("counts", [[[1, 2]]])
    def test_two_dimensional_counts(self, counts):
        _make_invalid(np_histogram=(counts, [0.0, 1.0, 2.0]))

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

    def test_nesting_with_num_bins_conflict(self):
        inner = Histogram([1, 2, 3], bins=2)
        with pytest.raises(ValueError, match="nesting"):
            Histogram(inner, num_bins=3)

    def test_nesting_other_transform_data_marks_invalid(self):
        h = _make_invalid(Text("hello"))  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="Invalid histogram input"):
            h.transform()


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
