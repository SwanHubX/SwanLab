"""
@author: mx
@file: eval.py
@time: 2026/11/25
@description: SwanLab 运行时评估 API
"""

from typing import Any, Dict, List, Optional

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.eval import InitEvalType

__all__ = ["upsert_eval"]


def upsert_eval(
    username: str,
    project: str,
    *,
    cuid: str,
    name: str,
    tags: Optional[List[str]] = None,
    link: Optional[str] = None,
    experiment_cuid: Optional[str] = None,
    step: Optional[int] = None,
    config: Optional[Dict[str, Any]] = None,
    retries: Optional[int] = None,
) -> InitEvalType:
    """
    创建或更新评估。

    指定 cuid 时服务端按 cuid 幂等：已存在则更新，不存在则以该 cuid 创建。
    config 传值时整体替换评估配置，省略时保留原配置。

    :param username: 所属用户名
    :param project: 所属项目名称
    :param cuid: 评估唯一标识符，由 SDK 生成
    :param name: 评估名称
    :param tags: 评估标签
    :param link: 外部链接
    :param experiment_cuid: 关联实验的 cuid 或 slug
    :param step: 被评 checkpoint 的训练步号
    :param config: 评估配置，格式为 {key: {value, desc, sort}}
    :param retries: 覆盖会话默认重试次数；None 使用会话默认值
    :return: 服务端保存的完整评估对象
    """
    body: Dict[str, Any] = {"cuid": cuid, "name": name}
    if tags:
        body["tags"] = list(tags)
    if link:
        body["link"] = link
    if experiment_cuid is not None:
        body["experimentCuid"] = experiment_cuid
    if step is not None:
        # 后端把 step 视为可为空的十进制字符串
        body["step"] = str(step)
    if config is not None:
        body["config"] = config
    resp = client.post(f"/evals/{username}/{project}/upsert", body, retries=retries)
    return resp.data
