"""
@author: cunyue
@file: batch.py
@time: 2025/12/4 13:50
@description: 分批上传函数
"""

import time
from typing import Union, Literal, List, TypedDict, Dict

from swanlab.core_python import get_client
from swanlab.log import swanlog


# 上传指标数据
class MetricDict(TypedDict):
    projectId: str
    experimentId: str
    type: str
    metrics: List[dict]
    flagId: Union[str, None]


def create_data(metrics: List[dict], metrics_type: str) -> MetricDict:
    """
    携带上传日志的指标信息
    """
    client = get_client()
    # Move 等实验需要将数据上传到根实验上
    exp_id = client.exp.root_exp_cuid or client.exp.cuid
    proj_id = client.exp.root_proj_cuid or client.proj.cuid
    assert proj_id is not None, "Project ID is empty."
    assert exp_id is not None, "Experiment ID is empty."
    flag_id = client.exp.flag_id
    return {
        "projectId": proj_id,
        "experimentId": exp_id,
        "type": metrics_type,
        "metrics": metrics,
        "flagId": flag_id,
    }


def _generate_chunks(data: Union[MetricDict, Dict, List], per_request_len: int):
    """
    生成器：统一处理字典和列表的分片逻辑
    yield: 分片后的数据块
    """
    # 情况1: 不分批
    if per_request_len == -1:
        yield data
        return

    # 情况2: 字典分批
    if isinstance(data, dict):
        metrics = data.get('metrics', [])
        # 如果 metrics 为空或长度不足，视为不需要分片的一整块
        if len(metrics) <= per_request_len:
            yield data
        else:
            for i in range(0, len(metrics), per_request_len):
                yield {
                    **data,
                    "metrics": metrics[i : i + per_request_len],
                }

    # 情况3: 列表分批
    elif isinstance(data, list):
        if len(data) <= per_request_len:
            yield data
        else:
            for i in range(0, len(data), per_request_len):
                yield data[i : i + per_request_len]


def trace_metrics(
    url: str,
    data: Union[MetricDict, list] = None,
    method: Literal['post', 'put'] = 'post',
    per_request_len: int = 1000,
):
    """
    分片指标上传方法
    """
    # 判断是否开启了分片模式（用于决定是否 sleep）
    # 这里的逻辑是：如果 per_request_len 不是 -1，且数据量确实超过了限制，则认为是分片模式
    is_split_mode = False
    if per_request_len != -1:
        total_len = len(data.get('metrics', [])) if isinstance(data, dict) else len(data or [])
        if total_len == 0:
            return
        is_split_mode = total_len > per_request_len
    client = get_client()
    # 遍历生成器产生的每一个数据块
    for chunk in _generate_chunks(data, per_request_len):
        # TODO: 暂时注释掉前置检查
        # 如果在发送过程中 client 变成了 pending，则中断后续发送
        # if client.pending:
        #     break

        # 调用被装饰的发送函数
        _, resp = getattr(client, method)(url, chunk)
        # 后置检查
        if resp and resp.status_code == 202:
            client.pending = True
            swanlog.warning(f"Client set to pending due to 202 response: {url}")
        # 分批发送时需要 sleep
        is_split_mode and time.sleep(1)
