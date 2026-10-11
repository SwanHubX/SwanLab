from typing import List, Literal, Sequence, Tuple, cast

import requests
from typing_extensions import deprecated

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.core_python.transport.sender.utils import transport_failure
from swanlab.sdk.internal.pkg import adapter
from swanlab.sdk.typings.core_python.api.upload import DeprecatedUploadColumn, DeprecatedUploadColumns
from swanlab.sdk.typings.core_python.transport import UploadResult


@deprecated("legacy projects without views will no longer be supported in v0.11.")
class DeprecatedHttpRecordSender(HttpRecordSender):
    """
    DeprecatedHttpRecordSender 是一个过时的类，用于适应多视图版本前的列上传逻辑
    """

    def upload_column(self, records: Sequence[Record], batch_size: int = 3000) -> UploadResult:
        encoded: List[Tuple[Record, DeprecatedUploadColumn]] = [
            (record, encode_desprecated_column(record.column)) for record in records
        ]
        confirmed: List[Record] = []
        for start in range(0, len(encoded), batch_size):
            batch = encoded[start : start + batch_size]
            try:
                upload_deprecated_columns(self._experiment_id, columns=[column for _, column in batch])
            except requests.RequestException as e:
                failure = transport_failure(
                    [record for record, _ in batch], e, source="swanlab_api", stage="submit_metrics"
                )
                untried = [record for record, _ in encoded[start + len(batch) :]]
                return UploadResult(
                    confirmed_records=confirmed,
                    remaining_records=[record for record, _ in batch] + untried,
                    failures=[failure],
                )
            confirmed.extend(record for record, _ in batch)
        return UploadResult(confirmed_records=confirmed)


def upload_deprecated_columns(experiment_id: str, *, columns: DeprecatedUploadColumns) -> None:
    """
    上传列信息
    """
    client.post(f"/experiment/{experiment_id}/columns", columns, retries=0)


def encode_desprecated_column(record: ColumnRecord) -> DeprecatedUploadColumn:
    """
    将列记录编码为后端所需的格式（DTO）
    NOTE: 此接口为多视图之前的版本
    """
    column: DeprecatedUploadColumn = {"key": record.column_key, "type": adapter.column[record.column_type]}  # noqa: F821
    # class: 是否为系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["class"] = cast(Literal["SYSTEM"], "SYSTEM")
    # name: 列的展示名称
    if record.column_name:
        column["name"] = record.column_name
    # section name: section 名称
    if record.section_name:
        column["sectionName"] = record.section_name
    # section type: section 类型，目前仅处理系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["sectionType"] = cast(Literal["SYSTEM"], "SYSTEM")
    # yRange: 数值列的 y 轴范围
    if record.column_type == ColumnType.COLUMN_TYPE_SCALAR:
        if record.HasField("y_range"):
            y_range = record.y_range
            min_value = y_range.min if y_range.HasField("min") else None
            max_value = y_range.max if y_range.HasField("max") else None
            column["yRange"] = (min_value, max_value)
    # chartName: 图表名称
    if record.chart_name:
        column["chartName"] = record.chart_name
    # chartIndex: 图表索引
    if record.chart_index:
        column["chartIndex"] = record.chart_index
    # metricName: 指标名称
    if record.metric_name:
        column["metricName"] = record.metric_name
    # metricColors: 指标颜色
    if record.HasField("metric_colors"):
        column["metricColors"] = (record.metric_colors.light, record.metric_colors.dark)
    return column


__all__ = ["DeprecatedHttpRecordSender"]
