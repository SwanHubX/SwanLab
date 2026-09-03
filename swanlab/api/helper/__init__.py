"""
@author: caddiesnew
@file: __init__.py
@time: 2026/9/3
@description: SwanLab OpenAPI 实体层公共辅助函数集合（统一出口，详见各领域模块）
"""

from swanlab.api.helper.extractor import (
    RELATIVE_TIME_AXIS,
    SCALAR_STATISTIC_FIELDS,
    STEP_AXIS,
    TIME_AXIS,
    align_entries_by_key,
    axis_request_params,
    builtin_x_axis,
    extract_first,
    extract_value_stats,
    merge_value_stats,
    stream_export_csv,
)
from swanlab.api.helper.request import (
    build_column_ref,
    build_export_payload,
    build_media_items,
    build_media_payload,
    build_scalar_payload,
    fetch_file_presigned_urls,
    fetch_presigned_urls,
)
from swanlab.api.helper.utils import (
    get_properties,
    parse_column_data_type,
    parse_timestamp_ms,
    resolve_run_path,
    strip_dict,
    validate_api_path,
    validate_column_params,
    validate_filter,
    validate_group,
    validate_metric_keys,
    validate_metric_log_level,
    validate_metric_type,
    validate_non_empty_string,
    validate_project_name,
    validate_sort,
    validate_update_active,
    validate_visibility,
    validate_x_axis,
)

__all__ = [
    # extractor —— 标量数据提取与解析（轴解析、统计合并、CSV 流式解析）
    "SCALAR_STATISTIC_FIELDS",
    "RELATIVE_TIME_AXIS",
    "STEP_AXIS",
    "TIME_AXIS",
    "align_entries_by_key",
    "axis_request_params",
    "builtin_x_axis",
    "extract_first",
    "extract_value_stats",
    "merge_value_stats",
    "stream_export_csv",
    # request —— House 请求体构建与预签名资源获取
    "build_column_ref",
    "build_export_payload",
    "build_media_items",
    "build_media_payload",
    "build_scalar_payload",
    "fetch_file_presigned_urls",
    "fetch_presigned_urls",
    # utils —— 参数校验与通用工具
    "get_properties",
    "parse_column_data_type",
    "parse_timestamp_ms",
    "resolve_run_path",
    "strip_dict",
    "validate_api_path",
    "validate_column_params",
    "validate_filter",
    "validate_group",
    "validate_metric_keys",
    "validate_metric_log_level",
    "validate_metric_type",
    "validate_non_empty_string",
    "validate_project_name",
    "validate_sort",
    "validate_update_active",
    "validate_visibility",
    "validate_x_axis",
]
