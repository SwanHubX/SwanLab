"""
@author: caddiesnew
@file: extractor.py
@time: 2026/9/3
@description: 标量指标数据提取辅助 — X 轴解析/分组、统计值合并、导出 CSV 流式解析
"""

import math
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from swanlab.api.typings.common import ApiResponseType, RangeQuery
from swanlab.api.typings.metric import ApiMetricXAxisType
from swanlab.sdk.internal.pkg import console, safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client

# ---------------------------------------------------------------------------
# X 轴解析与分组
# ---------------------------------------------------------------------------

# 回显形式：描述 metrics[].index 的实际语义
STEP_AXIS: ApiMetricXAxisType = {"type": "step"}
TIME_AXIS: ApiMetricXAxisType = {"type": "SYSTEM", "key": "time"}
RELATIVE_TIME_AXIS: ApiMetricXAxisType = {"type": "SYSTEM", "key": "relative_time"}


def builtin_x_axis(x_axis: str) -> Optional[ApiMetricXAxisType]:
    """显式内置轴字面量 → 回显形式；自定义 key 字符串返回 None。"""
    return {
        "step": STEP_AXIS,
        "time": TIME_AXIS,
        "relative_time": RELATIVE_TIME_AXIS,
    }.get(x_axis)


def axis_request_params(axis: ApiMetricXAxisType) -> Tuple[str, str]:
    """解析后的轴 → House 折线请求参数 ``(xType, xKey)``。

    CUSTOM 轴走 ``xKey``（服务端 xKey != "" 时自动切换 range/join 语义，xType 不参与）；
    SYSTEM 内置轴走 ``xType``；step 轴即现状默认。
    """
    axis_type = axis.get("type", "step")
    if axis_type == "CUSTOM":
        return ("step", axis.get("key", ""))
    key = axis.get("key", "")
    if axis_type == "SYSTEM" and key in ("time", "relative_time"):
        return (key, "")
    return ("step", "")


def group_keys_by_axis(axes: Dict[str, ApiMetricXAxisType]) -> List[Tuple[ApiMetricXAxisType, List[str]]]:
    """按解析后的轴分组（保持 key 首次出现顺序）。"""
    groups: Dict[Tuple[str, str], Tuple[ApiMetricXAxisType, List[str]]] = {}
    for key, axis in axes.items():
        group_id = (axis.get("type", "step"), axis.get("key", ""))
        if group_id not in groups:
            groups[group_id] = (axis, [])
        groups[group_id][1].append(key)
    return list(groups.values())


# ---------------------------------------------------------------------------
# 响应对齐与统计值合并
# ---------------------------------------------------------------------------

SCALAR_STATISTIC_FIELDS = ("min", "max", "avg", "median", "latest")


def extract_first(resp: ApiResponseType) -> Optional[Dict[str, Any]]:
    """从列表型 API 响应中提取第一个元素，失败返回 None。"""
    if resp.ok and isinstance(resp.data, list) and resp.data:
        return resp.data[0]
    return None


def align_entries_by_key(entries: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """将后端返回的列表按 ``key`` 字段映射为 dict，应对后端可能省略列或乱序返回。"""
    indexed: Dict[str, Dict[str, Any]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        key = entry.get("key")
        if isinstance(key, str) and key:
            indexed[key] = entry
    return indexed


def merge_value_stats(
    step_list: List[Dict[str, Any]],
    time_list: List[Dict[str, Any]],
    keys: List[str],
) -> List[Dict[str, Any]]:
    """合并 step/time 两种 x_type 的 value stat 响应为 per-key 统计字典。"""
    step_by_key = align_entries_by_key(step_list)
    time_by_key = align_entries_by_key(time_list)
    merged: List[Dict[str, Any]] = []
    for key in keys:
        step_entry = step_by_key.get(key, {})
        time_entry = time_by_key.get(key, {})
        entry: Dict[str, Any] = {}
        for field in SCALAR_STATISTIC_FIELDS:
            step_val = step_entry.get(field)
            time_val = time_entry.get(field)
            if isinstance(step_val, dict):
                stat = dict(step_val)
                if isinstance(time_val, dict) and time_val.get("index") is not None:
                    stat["timestamp"] = time_val["index"]
                entry[field] = stat
            elif isinstance(time_val, dict):
                entry[field] = dict(time_val)
        merged.append(entry)
    return merged


def extract_value_stats(
    step_resp: ApiResponseType,
    time_resp: ApiResponseType,
    keys: List[str],
) -> List[Dict[str, Any]]:
    """从并发的 step/time 响应中提取并合并 value stats。"""
    step_list = step_resp.data if step_resp.ok and isinstance(step_resp.data, list) else []
    time_list = time_resp.data if time_resp.ok and isinstance(time_resp.data, list) else []
    return merge_value_stats(step_list, time_list, keys)


# ---------------------------------------------------------------------------
# 导出 CSV 流式解析
# ---------------------------------------------------------------------------


@safe.decorator(message="Failed to download CSV")
def stream_export_csv(
    client: "Client",
    url: str,
    keys: List[str],
    rq: Optional[RangeQuery] = None,
    timeout: int = 30,
    x_key: str = "",
) -> Optional[Dict[str, List[Dict[str, Any]]]]:
    """Stream-download wide-format export CSV and parse per-key rows.

    CSV layout (from ``POST /house/metrics/scalar/export``)::

        step, {exp}-{key1}_step, {exp}-{key1}_timestamp,
              {exp}-{key2}_step, {exp}-{key2}_timestamp, … [, {exp}-{x}_step, {exp}-{x}_timestamp]

    One row per step; columns are interleaved as ``(value, timestamp)`` pairs.
    The ``{exp}-{key}_step`` column actually holds the metric **value** despite
    its name — the suffix is a House naming convention.

    When ``x_key`` is given (custom x axis), it is exported as one extra trailing
    column. Because House pivots by ``GROUP BY step``, the x value sits on the
    same row as every y value — alignment is applied in this single streaming
    pass.

    **NaN 占位语义**（区别于采样路径）：pivot 行集是所有请求列 step 的并集，每个 key
    在每一行都产出一个点——单元格缺失、NaN、无法解析时以 ``float("nan")`` 占位（自定义
    x 缺失时 ``index`` 同样占位 NaN），保证多 key / x 轴之间按下标对齐；±Inf 照常保留。
    仅时间戳类过滤（``last`` / ``type="timestamp"``）会丢弃缺失 timestamp 的行（无 ts
    可判）。采样路径（House 服务端）仍是 INNER JOIN + ``NOT isNaN`` + LTTB，会丢弃
    NaN 点——两条路径的行集差异为已知且有意为之。
    """
    import csv
    import time
    from collections import deque

    resp = client._session.get(url, stream=True, timeout=timeout)
    resp.raise_for_status()
    resp.encoding = "utf-8"
    lines = resp.iter_lines(decode_unicode=True)
    next(lines, None)  # skip header — column order is known from ``keys`` (+ optional trailing x)

    n_keys = len(keys)
    # x 列追加在所有 y 列之后，占最后一个 (value, timestamp) 槽位
    x_value_col = 1 + n_keys * 2 if x_key else None
    tail_limit = rq.tail if rq is not None and rq.tail is not None else None
    rows_per_key: List[Any] = [deque(maxlen=tail_limit) if tail_limit is not None else [] for _ in range(n_keys)]

    last_start_ts: Optional[int] = None
    if rq is not None and rq.last is not None:
        last_start_ts = int(time.time() * 1000) - rq.last

    # CSV is ``ORDER BY step`` — safe to break once step exceeds the range end.
    # ``type="custom"`` filters on the (possibly non-monotonic) x value domain instead,
    # so the step-order early break must be disabled there; head's count-based early
    # stop below remains valid in every mode.
    step_end_bound: Optional[int] = None
    if rq is not None and rq.type not in ("timestamp", "custom") and last_start_ts is None and rq.end is not None:
        # step/timestamp 语义下校验器已保证非负整值（float 存储，int() 截断安全）
        step_end_bound = int(rq.end)

    head_limit = rq.head if rq is not None and rq.head is not None else None
    _warned_missing_ts = False

    for row in csv.reader(lines):
        if not row:
            continue
        try:
            step = int(row[0])
        except (ValueError, IndexError):
            continue

        if step_end_bound is not None and step > step_end_bound:
            break

        # --- 自定义 x 列：缺失/NaN → index 以 NaN 占位（不丢行，保持多 key 按下标对齐） ---
        x_value: Optional[float] = None
        if x_value_col is not None:
            raw_x = row[x_value_col] if x_value_col < len(row) else ""
            if raw_x:
                try:
                    parsed_x = float(raw_x)
                except ValueError:
                    parsed_x = None
                if parsed_x is not None and not math.isnan(parsed_x):
                    x_value = parsed_x

        # --- type="custom"：按 x 值域纯谓词过滤（不假设单调，无提前终止依据） ---
        # x 占位 NaN 的行不参与谓词（IEEE 比较恒 False 会使 NaN 行意外通过有界过滤，
        # 这里显式放行占位行）
        if rq is not None and rq.type == "custom" and x_value is not None:
            if rq.start is not None and x_value < rq.start:
                continue
            if rq.end is not None and x_value > rq.end:
                continue

        for i in range(n_keys):
            vc = 1 + i * 2  # value column index
            tc = 2 + i * 2  # timestamp column index
            if vc >= len(row):
                continue
            raw_val = row[vc]
            value: float
            try:
                # 缺失（空串）/ NaN / 无法解析 → float("nan") 占位；±Inf 照常保留
                value = float(raw_val) if raw_val else float("nan")
            except ValueError:
                value = float("nan")

            ts: Optional[int] = None
            if tc < len(row) and row[tc]:
                try:
                    ts = int(row[tc])
                except ValueError:
                    pass

            if rq is not None:
                # --- ``last`` mode: filter by timestamp >= (now - last) ---
                if last_start_ts is not None:
                    if ts is None:
                        if not _warned_missing_ts:
                            console.warning("CSV row missing `timestamp` column.")
                            _warned_missing_ts = True
                        continue
                    if ts < last_start_ts:
                        continue
                # --- timestamp range mode ---
                elif rq.type == "timestamp":
                    if ts is None:
                        if not _warned_missing_ts:
                            console.warning("CSV row missing `timestamp` column.")
                            _warned_missing_ts = True
                        continue
                    if rq.start is not None and ts < rq.start:
                        continue
                    if rq.end is not None and ts > rq.end:
                        continue
                # --- step range mode (end handled by step_end_bound break) ---
                # 仅 type="step"：custom 已按 x 值域过滤，不得再用 start 当 step 下界
                elif rq.type == "step":
                    if rq.start is not None and step < rq.start:
                        continue

            item: Dict[str, Any] = {"step": step, "value": value}
            if x_value_col is not None:
                # x 缺失/NaN 时以 NaN 占位，保持 index 槽位存在
                item["index"] = x_value if x_value is not None else float("nan")
            if ts is not None:
                item["timestamp"] = ts
            rows_per_key[i].append(item)

        # head early-stop: all keys collected enough rows
        if head_limit is not None and all(len(r) >= head_limit for r in rows_per_key):
            break

    result: Dict[str, List[Dict[str, Any]]] = {}
    for i, key in enumerate(keys):
        rows = list(rows_per_key[i])
        if head_limit is not None:
            rows = rows[:head_limit]
        result[key] = rows
    return result
