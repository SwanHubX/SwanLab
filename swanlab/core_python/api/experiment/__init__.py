"""
@author: cunyue
@file: experiment.py
@time: 2025/12/11 18:37
@description: 定义实验相关的后端API接口
"""

from typing import TYPE_CHECKING, Dict, List, Literal, Optional, Union

from swanlab.core_python.api.type import RunType

from .utils import (
    extract_file_payloads,
    parse_column_type,
    to_camel_case,
    unwrap_api_payload,
)

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


MULTIPART_THRESHOLD: int = 100 * 1024 * 1024
PART_SIZE = 10 * 1024 * 1024


def send_experiment_heartbeat(
    client: "Client",
    *,
    cuid: str,
    flag_id: str,
):
    """
    发送实验心跳，保持实验处于活跃状态
    :param client: 已登录的客户端实例
    :param cuid: 实验唯一标识符
    :param flag_id: 实验标记ID
    """
    client.post(f"/house/experiments/{cuid}/heartbeat", {"flagId": flag_id})


def update_experiment_state(
    client: "Client",
    *,
    username: str,
    projname: str,
    cuid: str,
    state: Literal['FINISHED', 'CRASHED', 'ABORTED'],
    finished_at: str = None,
):
    """
    更新实验状态，注意此接口会将客户端标记为 pending 状态，表示实验已结束
    :param client: 已登录的客户端实例
    :param username: 实验所属用户名
    :param projname: 实验所属项目名称
    :param cuid: 实验唯一标识符
    :param state: 实验状态
    :param finished_at: 实验结束时间，格式为 ISO 8601，如果不提供则使用当前时间
    """
    put_data = {
        "state": state,
        "finishedAt": finished_at,
        "from": "sdk",
    }
    put_data = {k: v for k, v in put_data.items() if v is not None}  # 移除值为None的键
    client.put(f"/project/{username}/{projname}/runs/{cuid}/state", put_data)
    client.pending = True


def get_project_experiments(
    client: "Client",
    *,
    path: str,
    filters: Dict[str, object] = None,
) -> Union[List[RunType], Dict[str, List[RunType]]]:
    """
    获取指定项目下的所有实验信息
    若有实验分组，则返回一个字典，使用时需递归展平实验数据
    :param client: 已登录的客户端实例
    :param path: 项目路径 username/project
    :param filters: 筛选实验的条件，可选。支持以下特殊 key：
        - 'group': 按分组名称筛选，值为字符串
        - 'tags': 按标签筛选，值为字符串列表
        - 'name': 按实验名筛选，值为字符串
        - 'username': 按创建人筛选，值为字符串
        - 'job_type': 按任务类型筛选，值为字符串
    """
    # 特殊筛选条件配置：用户侧 key -> 后端 key 和操作符
    special_filter_config = {
        "group": {"key": "cluster", "op": "EQ"},
        "tags": {"key": "labels", "op": "IN"},
        "name": {"key": "name", "op": "EQ"},
        "username": {"key": "user.username", "op": "EQ"},
        "job_type": {"key": "job", "op": "EQ"},
    }

    parsed_filters = []

    if filters:
        for key, value in filters.items():
            if key in special_filter_config:
                # 特殊字段处理
                config = special_filter_config[key]
                # tags 需要转换为列表
                filter_value = list(value) if key == "tags" and isinstance(value, (list, tuple)) else [value]
                parsed_filters.append(
                    {
                        "key": config["key"],
                        "active": True,
                        "value": filter_value,
                        "op": config["op"],
                        "type": 'STABLE',
                    }
                )
            else:
                # 常规字段处理
                parsed_filters.append(
                    {
                        "key": to_camel_case(key) if parse_column_type(key) == 'STABLE' else key.split('.', 1)[-1],
                        "active": True,
                        "value": [value],
                        "op": 'EQ',
                        "type": parse_column_type(key),
                    }
                )

    res = client.post(f"/project/{path}/runs/shows", data={'filters': parsed_filters})
    return res[0]


def get_single_experiment(client: "Client", *, path: str) -> RunType:
    """
    获取指定项目下的所有实验信息
    若有实验分组，则返回一个字典，使用时需递归展平实验数据
    :param client: 已登录的客户端实例
    :param path: 实验路径 username/project/expid
    """
    proj_path, expid = path.rsplit('/', 1)
    res = client.get(f"/project/{proj_path}/runs/{expid}")
    return res[0]


def get_experiment_metrics(client: "Client", *, expid: str, key: str) -> Dict[str, str]:
    """
    获取指定字段的指标数据，返回csv网址
    :param client: 已登录的客户端实例
    :param expid: 实验cuid
    :param key: 指定字段列表
    """
    res = client.get(f"/experiment/{expid}/column/csv", params={'key': key})
    return res[0]


def delete_experiment(client: "Client", *, path: str):
    """
    删除指定实验
    :param client: 已登录的客户端实例
    :param path: 实验路径 'username/project/expid'
    """
    proj_path, expid = path.rsplit('/', 1)
    client.delete(f"/project/{proj_path}/runs/{expid}")


def prepare_upload(client: "Client", exp_id: str, files: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """
    创建普通文件上传任务，返回预签名上传地址列表。
    """
    if len(files) == 0:
        return []
    data, _ = client.post(f"/experiment/{exp_id}/files/prepare", {"files": files})
    return extract_file_payloads(unwrap_api_payload(data))


def complete_upload(client: "Client", exp_id: str, names: List[str], state: str = "UPLOADED") -> None:
    """
    标记普通文件上传完成。
    """
    if len(names) == 0:
        return
    client.post(
        f"/experiment/{exp_id}/files/complete",
        {"files": [{"name": name, "state": state} for name in names]},
    )


def prepare_multipart(
    client: "Client",
    exp_id: str,
    name: str,
    size: int,
    part_count: int,
    md5: str,
    mime_type: Optional[str] = None,
) -> Dict[str, object]:
    """
    创建分片上传任务，返回上传地址和上传上下文。
    """
    file_payload: Dict[str, object] = {
        "name": name,
        "size": size,
        "md5": md5,
        "count": part_count,
    }
    if mime_type is not None:
        file_payload["mimeType"] = mime_type
    data, _ = client.post(
        f"/experiment/{exp_id}/files/prepare-multipart",
        {"files": [file_payload]},
    )
    payloads = extract_file_payloads(unwrap_api_payload(data))
    if len(payloads) == 0:
        raise ValueError("Multipart prepare API returned empty file payloads.")
    return payloads[0]


def complete_multipart(
    client: "Client",
    exp_id: str,
    name: str,
    upload_id: str,
    state: str = "UPLOADED",
) -> None:
    """
    标记分片上传完成，并通知后端执行合并。
    """
    client.post(
        f"/experiment/{exp_id}/files/complete-multipart",
        {"files": [{"name": name, "uploadId": upload_id, "state": state}]},
    )


__all__ = [
    "MULTIPART_THRESHOLD",
    "PART_SIZE",
    "send_experiment_heartbeat",
    "update_experiment_state",
    "get_project_experiments",
    "get_single_experiment",
    "get_experiment_metrics",
    "delete_experiment",
    "prepare_upload",
    "complete_upload",
    "prepare_multipart",
    "complete_multipart",
]
