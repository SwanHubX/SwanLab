"""
@author: Zirui Cai
@file: checkpoint.py
@time: 2026/11/25
@description: Checkpoint 登记 API
"""

from typing import List, Optional

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper


def register_checkpoint(
    *,
    experiment_id: str,
    uri: str,
    step: Optional[int] = None,
    name: Optional[str] = None,
    description: Optional[str] = None,
    metadata: Optional[dict] = None,
    labels: Optional[List[str]] = None,
) -> str:
    """登记一个 checkpoint，返回平台分配的 cuid。

    :param experiment_id: 来源实验的 cuid（即 run_id）
    :param uri: checkpoint 路径或 URI（本地路径、S3、HF Hub URL 等）
    :param step: 保存时的训练步数；None 表示未记录
    :param name: checkpoint 名称；None 时由 Server 从 URI 末段生成
    :param description: 备注说明
    :param metadata: 用户自定义元数据
    :param labels: 标签名列表，如 ["best", "ready-for-eval"]
    :return: 平台分配的 checkpoint cuid
    """
    labels_body = [{"name": n} for n in labels] if labels else None
    body = helper.strip_none(
        {
            "experimentId": experiment_id,
            "uri": uri,
            "step": step,
            "name": name,
            "description": description,
            "metadata": metadata,
            "labels": labels_body,
        }
    )
    # retries=0：重试机制交给 transport 层的退避调度，与 save.py 保持一致
    resp = client.post("/checkpoint", body, retries=0)
    return resp.data["cuid"]
