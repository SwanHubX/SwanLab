"""
@author: Nexisato
@file: view.py
@time: 2026/9/3
@description: DEFAULT 视图 X 轴自动发现辅助模块（进程内缓存 + 优雅回落）
"""

from typing import Any, Dict, List, Optional, Tuple

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.metric import ApiMetricXAxisType
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor

_DEFAULT_VIEW_NAME = "Default"

# (api_host, proj_path) → title → axis。Api 与 CLI 均短生命周期，不做失效策略；
# 长寿命进程中视图配置的会话内变更不会被感知（已知取舍，文档说明）。
_XAXIS_CACHE: Dict[Tuple[str, str], Dict[str, ApiMetricXAxisType]] = {}
# 发现失败的告警去重键（成功后清除，允许下次失败再次提示）
_WARNED: set = set()


class _ViewEntity(BaseEntity):
    """仅为复用 BaseEntity 的 _get 安全请求通道（失败不抛异常，返回 ApiResponseType）。"""

    def json(self) -> Dict[str, Any]:
        return {}


def fetch_default_view_xaxis_map(ctx: ApiClientContext, proj_path: str) -> Dict[str, ApiMetricXAxisType]:
    """解析 DEFAULT 视图中每个 scalar key 的 X 轴配置（title → axis）。

    任何失败（非 200 / 超时 / 响应体畸形）只 warning 一次并返回空映射 ——
    指标查询绝不因发现失败而失败，调用方回落 step。
    """
    cache_key = (ctx.api_host, proj_path)
    cached = _XAXIS_CACHE.get(cache_key)
    if cached is not None:
        return cached
    discovered = _discover(ctx, proj_path)
    if discovered is None:
        if cache_key not in _WARNED:
            console.warning(
                f"Failed to resolve the DEFAULT view x-axis config for '{proj_path}'; "
                "metrics will fall back to the 'step' axis."
            )
            _WARNED.add(cache_key)
        return {}
    _XAXIS_CACHE[cache_key] = discovered
    _WARNED.discard(cache_key)
    return discovered


def _as_list(data: Any) -> Optional[List[Any]]:
    """容忍 list 或 ``{"list": [...]}`` 两种响应形态；畸形返回 None。"""
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and isinstance(data.get("list"), list):
        return data["list"]
    return None


def _discover(ctx: ApiClientContext, proj_path: str) -> Optional[Dict[str, ApiMetricXAxisType]]:
    """执行发现请求；失败返回 None（区别于"无视图"的合法空映射）。"""
    entity = _ViewEntity(ctx)

    views_resp = entity._get(f"/views/{proj_path}")
    views = _as_list(views_resp.data) if views_resp.ok else None
    if views is None:
        return None
    view_dicts = [v for v in views if isinstance(v, dict)]
    view = next((v for v in view_dicts if v.get("name") == _DEFAULT_VIEW_NAME), None)
    if view is None and view_dicts:
        # 兜底：第一个视图
        view = view_dicts[0]
    index = view.get("index") if view is not None else None
    if index is None:
        # 无视图 → 空映射（合法状态，缓存）
        return {}

    sections_resp = entity._get(f"/sections/{proj_path}/{index}")
    sections = _as_list(sections_resp.data) if sections_resp.ok else None
    if sections is None:
        return None
    section_indices = [s.get("index") for s in sections if isinstance(s, dict) and s.get("index") is not None]
    if not section_indices:
        return {}

    # 逐 section 请求 charts（含 config），并发拉取；Hidden section 中的图表同样计入
    with SafeThreadPoolExecutor(max_workers=min(len(section_indices), 4)) as pool:
        futures = [pool.submit(entity._get, f"/charts/{proj_path}/{index}/{si}") for si in section_indices]
        chart_resps = [f.result() for f in futures]

    axis_map: Dict[str, ApiMetricXAxisType] = {}
    for resp in chart_resps:
        charts = _as_list(resp.data) if resp.ok else None
        if charts is None:
            return None
        for chart in charts:
            if not isinstance(chart, dict):
                continue
            axis = _parse_chart_axis(chart)
            if axis is not None:
                axis_map[chart.get("title", "")] = axis
    return axis_map


def _parse_chart_axis(chart: Dict[str, Any]) -> Optional[ApiMetricXAxisType]:
    """从单张图表配置解析 X 轴。

    仅接受固定轴（``fixedXAxis`` 为真，服务端实际存放于 ``config`` 内，顶层值兼容）
    且 ``config.xAxis`` 存在的图表，且 ``config.yAxis`` 必须为单 key 且等于 title
    （DEFAULT 视图自动图满足 ``title === key``；用户拼接的多 key 图表跳过 —— 已记录
    的限制）。yAxis 实际形状为对象数组（取 ``key`` 字段），纯字符串数组 defensively 兼容。
    """
    title = chart.get("title")
    if not isinstance(title, str) or not title:
        return None
    config = chart.get("config")
    if not isinstance(config, dict):
        return None
    axis = config.get("xAxis")
    if not isinstance(axis, dict):
        return None
    if not (chart.get("fixedXAxis") or config.get("fixedXAxis")):
        return None
    y_axis = config.get("yAxis")
    if isinstance(y_axis, list):
        if len(y_axis) != 1:
            return None
        entry = y_axis[0]
        y_key = entry.get("key") if isinstance(entry, dict) else entry
        if y_key != title:
            return None
    elif y_axis != title:
        return None
    axis_class = axis.get("class") or axis.get("source") or ""
    key = axis.get("key") or ""
    if axis_class == "CUSTOM" and key:
        return {"type": "CUSTOM", "key": key}
    if axis_class == "SYSTEM":
        # DEFAULT 视图的 SYSTEM 轴即 relative_time
        return {"type": "SYSTEM", "key": key or "relative_time"}
    return None
