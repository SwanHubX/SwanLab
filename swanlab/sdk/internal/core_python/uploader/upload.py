"""
@author: caddiesnew
@file: upload.py
@time: 2026/3/19
@description: 定义各类型上传函数
"""

import json
from pathlib import Path
from typing import Callable, List, Optional, Union

import yaml

from swanlab.sdk.internal import adapter
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.system.v1.console_pb2 import StreamType
from swanlab.sdk.internal.pkg import console

from .batch import create_data, trace_metrics
from .model import FileModel

HOUSE_URL = "/house/metrics"
RecordLike = Union[Record, bytes]


def _load_record(record: RecordLike) -> Record:
    if isinstance(record, Record):
        return record
    if isinstance(record, bytes):
        parsed = Record()
        parsed.ParseFromString(record)
        return parsed
    raise TypeError(f"Unsupported record type: {type(record).__name__}")


def _timestamp_to_json(ts) -> Optional[str]:
    if ts is None or (ts.seconds == 0 and ts.nanos == 0):
        return None
    return ts.ToJsonString()


def _media_item_to_dict(item) -> dict:
    payload = {"filename": item.filename, "sha256": item.sha256}
    if hasattr(item, "size"):
        payload["size"] = item.size
    if hasattr(item, "caption") and item.caption != "":
        payload["caption"] = item.caption
    return payload


def _metric_record_to_dict(record_like: Union[dict, RecordLike]) -> dict:
    if isinstance(record_like, dict):
        return record_like

    record = _load_record(record_like)
    metric = record.metric
    value_type = metric.WhichOneof("value")
    if value_type is None:
        raise ValueError("Metric record has no value field.")

    payload = {
        "key": metric.key,
        "index": metric.step,
        "type": adapter.column_type.get(metric.type),
    }
    timestamp = _timestamp_to_json(metric.timestamp)
    if timestamp is not None:
        payload["timestamp"] = timestamp

    if value_type == "scalar":
        payload["scalar"] = {"number": metric.scalar.number}
        return payload

    items = [_media_item_to_dict(item) for item in getattr(metric, value_type).items]
    payload[value_type] = {"items": items}
    return payload


def _column_record_to_dict(column_like: Union[dict, RecordLike]) -> dict:
    if isinstance(column_like, dict):
        return column_like

    record = _load_record(column_like)
    column = record.column
    payload = {
        "class": column.column_class,
        "type": adapter.column_type.get(column.column_type),
        "key": column.column_key,
    }
    if column.column_name:
        payload["name"] = column.column_name
    if column.section_name:
        payload["sectionName"] = column.section_name
    if column.section_type:
        payload["sectionType"] = column.section_type
    if column.HasField("y_range"):
        payload["yRange"] = {"min": column.y_range.min, "max": column.y_range.max}
    if column.chart_name:
        payload["chartName"] = column.chart_name
    if column.chart_index:
        payload["chartIndex"] = column.chart_index
    if column.metric_name:
        payload["metricName"] = column.metric_name
    if len(column.metric_colors) > 0:
        payload["metricColors"] = list(column.metric_colors)
    return payload


def _log_record_to_dict(log_like: Union[dict, RecordLike]) -> dict:
    if isinstance(log_like, dict):
        return log_like

    record = _load_record(log_like)
    console_record = record.console
    payload = {
        "line": console_record.line,
        "stream": console_record.stream,
        "level": "ERROR" if console_record.stream == StreamType.STREAM_TYPE_STDERR else "INFO",
    }
    timestamp = _timestamp_to_json(console_record.timestamp)
    if timestamp is not None:
        payload["timestamp"] = timestamp
    return payload


def _read_text_file(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8")


def _read_json_file(path: Path) -> Optional[dict]:
    text = _read_text_file(path)
    if text is None:
        return None
    if text.strip() == "":
        return {}
    return json.loads(text)


def _read_yaml_file(path: Path) -> Optional[dict]:
    text = _read_text_file(path)
    if text is None:
        return None
    if text.strip() == "":
        return {}
    data = yaml.safe_load(text)
    return data if data is not None else {}


def _file_model_from_records(files: List[RecordLike], files_dir: Optional[Path]) -> FileModel:
    if files_dir is None:
        raise RuntimeError("files_dir is required when uploading file records.")

    file_types = {_load_record(file_record).WhichOneof("record_type") for file_record in files}
    return FileModel(
        requirements=_read_text_file(files_dir / adapter.filename.requirements) if "requirements" in file_types else None,
        metadata=_read_json_file(files_dir / adapter.filename.metadata) if "metadata" in file_types else None,
        config=_read_yaml_file(files_dir / adapter.filename.config) if "config" in file_types else None,
        conda=_read_text_file(files_dir / adapter.filename.conda) if "conda" in file_types else None,
    )


def upload_logs(logs: List[Union[dict, RecordLike]], upload_callback: Optional[Callable] = None) -> None:
    """
    上传日志信息（ConsoleRecord / 序列化后的 Record）。
    """
    if len(logs) == 0:
        return console.debug("No logs to upload.")
    metrics = [_log_record_to_dict(log) for log in logs]
    data = create_data(metrics, "log")
    trace_metrics(HOUSE_URL, data, upload_callback=upload_callback)


def upload_scalar_metrics(
    scalar_metrics: List[Union[dict, RecordLike]], upload_callback: Optional[Callable] = None
) -> None:
    """
    上传标量指标数据（DataRecord scalar / 序列化后的 Record）。
    """
    if len(scalar_metrics) == 0:
        return console.debug("No scalar metrics to upload.")
    data = create_data([_metric_record_to_dict(metric) for metric in scalar_metrics], "scalar")
    trace_metrics(HOUSE_URL, data, upload_callback=upload_callback)


def upload_media_metrics(
    media_metrics: List[Union[dict, RecordLike]], upload_callback: Optional[Callable] = None
) -> None:
    """
    上传媒体指标数据（DataRecord media / 序列化后的 Record）。
    """
    if len(media_metrics) == 0:
        return console.debug("No media metrics to upload.")
    # TODO: 媒体文件 buffer 上传（upload_to_cos）待后续接入
    data = create_data([_metric_record_to_dict(metric) for metric in media_metrics], "media")
    trace_metrics(HOUSE_URL, data, upload_callback=upload_callback)


def upload_columns(columns: List[Union[dict, RecordLike]], upload_callback: Optional[Callable] = None) -> None:
    """
    批量上传并创建 columns。
    """
    if len(columns) == 0:
        return console.debug("No columns to upload.")
    trace_metrics(
        "/experiment/columns",
        [_column_record_to_dict(column) for column in columns],
        per_request_len=3000,
        upload_callback=upload_callback,
    )


def upload_files(files: List[Union[FileModel, RecordLike]], files_dir: Optional[Path] = None) -> None:
    """
    上传 files 文件夹中的内容（metadata/requirements/conda/config 聚合）。
    """
    if len(files) == 0:
        return console.debug("No files to upload.")
    if all(isinstance(file_model, FileModel) for file_model in files):
        file_model = FileModel.create(files)
    else:
        file_model = _file_model_from_records(files, files_dir=files_dir)
    if file_model.empty:
        return
    data = file_model.to_dict()
    trace_metrics("/profile", data, method="put", per_request_len=-1)


# 上传类型 → 上传函数 + 是否支持回调
UPLOAD_DISPATCH = {
    "scalar": {"upload": upload_scalar_metrics, "has_callback": True},
    "media": {"upload": upload_media_metrics, "has_callback": True},
    "column": {"upload": upload_columns, "has_callback": True},
    "log": {"upload": upload_logs, "has_callback": True},
    "file": {"upload": upload_files, "has_callback": False},
}


__all__ = [
    "upload_logs",
    "upload_scalar_metrics",
    "upload_media_metrics",
    "upload_columns",
    "upload_files",
    "UPLOAD_DISPATCH",
]
