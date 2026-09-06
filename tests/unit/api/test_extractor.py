"""
@author: caddiesnew
@time: 2026/9/3
@description: swanlab/api/helper/extractor.py 核心流式 CSV 解析与对齐单测
"""

import math
from unittest.mock import MagicMock

from swanlab.api.helper.extractor import stream_export_csv
from swanlab.api.typings.common import RangeQuery


def _mock_client(csv_content: str):
    mock_resp = MagicMock()
    mock_resp.raise_for_status.return_value = None
    mock_resp.encoding = "utf-8"
    mock_resp.iter_lines.return_value = iter(csv_content.strip().splitlines())

    mock_client = MagicMock()
    mock_client._session.get.return_value = mock_resp
    return mock_client


class TestStreamExportCsv:
    def test_non_custom_x_axis_skips_missing_cells(self):
        """非自定义轴（纯 step 轴，not x_key）：
        旧契约守恒，空单元格或无效值直接跳过，不填充 NaN。
        """
        csv_text = (
            "step,c1,c1_timestamp,c2,c2_timestamp\n"
            "0,0.5,1700000000000,0.8,1700000000000\n"
            "1,0.4,1700000001000,,1700000001000\n"
            "2,,1700000002000,0.9,1700000002000\n"
        )
        client = _mock_client(csv_text)
        result = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["loss", "acc"],
            rq=None,
            x_key="",
        )
        assert result is not None
        loss_metrics = result["loss"]
        acc_metrics = result["acc"]

        assert len(loss_metrics) == 2
        assert [m["step"] for m in loss_metrics] == [0, 1]
        assert [m["value"] for m in loss_metrics] == [0.5, 0.4]
        assert not any(math.isnan(m["value"]) for m in loss_metrics)

        assert len(acc_metrics) == 2
        assert [m["step"] for m in acc_metrics] == [0, 2]
        assert [m["value"] for m in acc_metrics] == [0.8, 0.9]
        assert not any(math.isnan(m["value"]) for m in acc_metrics)

    def test_custom_x_axis_nan_placeholders_and_index(self):
        """自定义 x 轴（x_key 存在且不在 keys 中）：
        x_key 作为追加列位于末尾。缺失单元格以 NaN 占位，index 为对应的 x_value。
        """
        csv_text = (
            "step,c1,c1_ts,c2,c2_ts,c3,c3_ts\n"
            "0,0.5,1000,0.8,1000,1.0,1000\n"
            "1,0.4,2000,,2000,2.0,2000\n"
            "2,,3000,0.9,3000,3.0,3000\n"
        )
        client = _mock_client(csv_text)
        result = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["loss", "acc"],
            rq=None,
            x_key="epoch",
        )
        assert result is not None
        loss_metrics = result["loss"]
        acc_metrics = result["acc"]

        assert len(loss_metrics) == 3
        assert len(acc_metrics) == 3

        assert [m["index"] for m in loss_metrics] == [1.0, 2.0, 3.0]
        assert [m["index"] for m in acc_metrics] == [1.0, 2.0, 3.0]

        assert math.isnan(acc_metrics[1]["value"])
        assert math.isnan(loss_metrics[2]["value"])

    def test_x_key_in_keys_column_reuse_and_identity(self):
        """当 x_key 存在于 keys 中（如 keys=["epoch", "loss"], x_key="epoch"）：
        列复用，不追加末尾列，且 epoch 序列的 value == index。
        """
        csv_text = "step,c1,c1_ts,c2,c2_ts\n0,10.0,1000,0.5,1000\n1,20.0,2000,0.4,2000\n"
        client = _mock_client(csv_text)
        result = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["epoch", "loss"],
            rq=None,
            x_key="epoch",
        )
        assert result is not None
        epoch_metrics = result["epoch"]
        loss_metrics = result["loss"]

        for m in epoch_metrics:
            assert m["value"] == m["index"]
        assert [m["index"] for m in epoch_metrics] == [10.0, 20.0]

        assert [m["index"] for m in loss_metrics] == [10.0, 20.0]
        assert [m["value"] for m in loss_metrics] == [0.5, 0.4]

    def test_custom_range_bounds_and_nan_drop(self):
        """type="custom" 时：
        按 x_value 过滤 start/end，支持负值范围；
        当 x 缺失为 NaN/None 时，在有界查询下被过滤丢弃。
        """
        csv_text = (
            "step,c1,c1_ts,c2,c2_ts\n0,0.1,1000,-1.0,1000\n1,0.2,2000,,2000\n2,0.3,3000,0.5,3000\n3,0.4,4000,2.0,4000\n"
        )
        client = _mock_client(csv_text)
        rq = RangeQuery(type="custom", start=-0.5, end=1.0)
        result = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["loss"],
            rq=rq,
            x_key="lr",
        )
        assert result is not None
        metrics = result["loss"]
        assert len(metrics) == 1
        assert metrics[0]["index"] == 0.5
        assert metrics[0]["value"] == 0.3

    def test_range_head_and_tail_with_custom_x(self):
        """head 和 tail 在 custom x 轴下的截断行为。"""
        csv_text = (
            "step,c1,c1_ts,c2,c2_ts\n"
            "0,0.1,1000,1.0,1000\n"
            "1,0.2,2000,2.0,2000\n"
            "2,0.3,3000,3.0,3000\n"
            "3,0.4,4000,4.0,4000\n"
        )
        # 测试 head=2
        client = _mock_client(csv_text)
        result_head = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["loss"],
            rq=RangeQuery(type="custom", head=2),
            x_key="epoch",
        )
        assert result_head is not None
        assert [m["index"] for m in result_head["loss"]] == [1.0, 2.0]

        # 测试 tail=2
        client = _mock_client(csv_text)
        result_tail = stream_export_csv(
            client=client,
            url="https://mock/csv",
            keys=["loss"],
            rq=RangeQuery(type="custom", tail=2),
            x_key="epoch",
        )
        assert result_tail is not None
        assert [m["index"] for m in result_tail["loss"]] == [3.0, 4.0]
