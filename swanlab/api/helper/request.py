"""
@author: caddiesnew
@file: request.py
@time: 2026/9/3
@description: House 查询请求辅助 — 请求体构建、预签名链接获取与媒体项构建
"""

from typing import Any, Dict, List, Optional

from swanlab.api.base import BaseEntity
from swanlab.api.typings.metric import ApiMediaItemDataType

# ---------------------------------------------------------------------------
# 请求体构建（纯函数）
# ---------------------------------------------------------------------------


def build_column_ref(
    experiment_id: str,
    created_at: int,
    key: str,
    root_pro_id: str = "",
    root_exp_id: str = "",
) -> Dict[str, Any]:
    # createdAt 为 House 查询的数据入库时间下界，必传
    ref: Dict[str, Any] = {"experimentId": experiment_id, "key": key, "createdAt": created_at}
    if root_pro_id:
        ref["rootProId"] = root_pro_id
    if root_exp_id:
        ref["rootExpId"] = root_exp_id
    return ref


def build_scalar_payload(
    project_id: str,
    run_id: str,
    created_at: int,
    keys: List[str],
    sample: int = 1500,
    x_type: str = "step",
    x_key: Optional[str] = None,
    root_pro_id: str = "",
    root_exp_id: str = "",
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "projectId": project_id,
        "xType": x_type,
        "range": [0, 0],
        "columns": [build_column_ref(run_id, created_at, key, root_pro_id, root_exp_id) for key in keys],
        "num": sample if sample <= 1500 else 1500,
    }
    # xKey != "" 时服务端自动把 range/join 语义切换为 x 值域；"range": [0, 0] 用不到
    if x_key:
        payload["xKey"] = x_key
    return payload


def build_media_payload(
    project_id: str,
    run_id: str,
    created_at: int,
    keys: List[str],
    step: Optional[int] = None,
    root_pro_id: str = "",
    root_exp_id: str = "",
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "projectId": project_id,
        "columns": [build_column_ref(run_id, created_at, key, root_pro_id, root_exp_id) for key in keys],
    }
    if step is not None:
        payload["step"] = step
    return payload


def build_export_payload(
    project_id: str,
    run_id: str,
    created_at: int,
    keys: List[str],
    experiment_name: str = "",
    root_pro_id: str = "",
    root_exp_id: str = "",
) -> Dict[str, Any]:
    """Build payload for ``POST /house/metrics/scalar/export``.

    ``experimentName`` is required by the API but only used for CSV column
    headers — the actual query uses ``experimentId``. Falls back to ``run_id``
    when the real name is unavailable.
    """
    exp_name = experiment_name or run_id
    columns: List[Dict[str, Any]] = []
    for key in keys:
        # createdAt 为 House 查询的数据入库时间下界，必传
        col: Dict[str, Any] = {
            "experimentName": exp_name,
            "experimentId": run_id,
            "key": key,
            "createdAt": created_at,
        }
        if root_pro_id:
            col["rootProId"] = root_pro_id
        if root_exp_id:
            col["rootExpId"] = root_exp_id
        columns.append(col)
    return {"projectId": project_id, "columns": columns}


# ---------------------------------------------------------------------------
# 预签名链接获取与媒体项构建
# ---------------------------------------------------------------------------


def fetch_presigned_urls(entity: BaseEntity, prefix: str, paths: List[str]) -> Dict[str, str]:
    """批量获取预签名下载链接，返回 path → url 映射。"""
    if not paths:
        return {}
    resp = entity._post("/resources/presigned/get", data={"prefix": prefix, "paths": paths})
    if not resp.ok or not isinstance(resp.data, dict):
        return {}
    urls = resp.data.get("urls") or []
    return dict(zip(paths, urls)) if urls else {}


def fetch_file_presigned_urls(entity: BaseEntity, paths: List[str]) -> Dict[str, str]:
    """通过完整资源路径批量获取预签名下载链接，返回 path → url 映射。"""
    if not paths:
        return {}
    resp = entity._post("/files/presigned/get", data={"paths": paths})
    if not resp.ok or not isinstance(resp.data, dict):
        return {}
    urls = resp.data.get("urls") or []
    return dict(zip(paths, urls)) if urls else {}


def build_media_items(
    entry: Dict[str, Any],
    url_map: Dict[str, str],
) -> List[ApiMediaItemDataType]:
    """将单个 metric entry 的 data/more 合并为 items，注入预签名 url。"""
    # 后端对"无数据"可能返回显式 null（dict.get 的默认值不生效），统一兜底为空列表
    paths = entry.get("data") or []
    mores = entry.get("more") or []
    items: List[ApiMediaItemDataType] = []
    for i, path in enumerate(paths):
        item: ApiMediaItemDataType = {}
        if path in url_map:
            item["url"] = url_map[path]
        if i < len(mores) and isinstance(mores[i], dict):
            item.update(mores[i])
        items.append(item)
    return items
