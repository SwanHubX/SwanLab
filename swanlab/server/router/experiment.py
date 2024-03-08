from fastapi import APIRouter
from ..module.resp import PARAMS_ERROR_422

from ..controller.experiment import (
    # 获取实验信息
    get_experiment_info,
    # 获取 tag 信息
    get_tag_data,
    # 获取实验状态
    get_experiment_status,
    # 获取实验的总结数据
    get_experiment_summary,
    # 获取实验最近日志
    get_recent_logs,
    # 获取实验图标
    get_experimet_charts,
    # 更改实验信息
    update_experiment_info,
    # 删除实验
    delete_experiment,
    # 停止实验
    stop_experiment,
    # 获取实验依赖
    get_experiment_requirements,
    # 修改实验可见性
    change_experiment_visibility,
)
from fastapi import Request
from urllib.parse import quote

router = APIRouter()


@router.get("/{experiment_id}")
def _(experiment_id: int):
    """获取实验信息

    Parameters
    ----------
    experiment_id : int
        实验唯一 id

    Returns
    -------
    dict
        实验信息
    """

    return get_experiment_info(experiment_id)


# COMPAT 由于fastapi不支持%2F的路径转换，所以采用通配符:path，并且在下面将path进行quote编码
@router.get("/{experiment_id}/tag/{tag:path}")
def _(experiment_id: int, tag: str) -> dict:
    """获取表单数据

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，使用时需要 URIComponent 编码
    """
    tag = quote(tag, safe="")

    return get_tag_data(experiment_id, tag)


@router.get("/{experiment_id}/status")
def _(experiment_id: int):
    """获取实验状态以及实验图表配置，用于实时更新实验状态

    Parameters
    ----------
    experiment_id : int
        实验唯一id，路径传参
    """

    return get_experiment_status(experiment_id)


@router.get("/{experiment_id}/summary")
def _(experiment_id: int) -> dict:
    """获取实验的总结数据——每个tag的最后一个setp的data

    Parameters
    ----------
    experiment_id : int
        实验id

    Returns
    -------
    dict
        响应信息
    """

    return get_experiment_summary(experiment_id)


@router.get("/{experiment_id}/recent_log")
async def _(experiment_id: int):
    """一下返回最多 MAX_NUM 条打印记录

    Parameters
    ----------
    experiment_id : int
        实验唯一ID
    """

    return get_recent_logs(experiment_id)


@router.get("/{experiment_id}/chart")
def _(experiment_id: int):
    """获取图标信息

    Parameters
    ----------
    experiment_id : int
        实验唯一 ID
    """

    return get_experimet_charts(experiment_id)


@router.patch("/{experiment_id}")
async def _(experiment_id: int, request: Request):
    """修改实验的信息

    Parameters
    ----------
    experiment_id : int
        实验id
    body : Body
        name: str
            实验名称
        description: str
            实验描述

    Returns
    -------
    dict
    """

    return await update_experiment_info(experiment_id, request)


@router.delete("/{experiment_id}")
def _(experiment_id: int):
    """删除实验

    注意，需要先判断当前实验是否正在运行中，不可删除运行中的实验

    Parameters
    ----------
    experiment_id : Int
        实验唯一ID

    Returns
    -------
    project : Dictionary
        删除实验后的项目信息，提供给前端更新界面
    """

    return delete_experiment(experiment_id)


@router.get("/{experiment_id}/stop")
def _(experiment_id: int):
    """停止实验

    Parameters
    ----------
    experiment_id : Int
        实验唯一ID
    """

    return stop_experiment(experiment_id)


# 获取实验依赖
@router.get("/{experiment_id}/requirements")
def _(experiment_id: int):
    """获取实验依赖

    Parameters
    ----------
    experiment_id : int
        实验唯一ID

    Returns
    -------
    requirements: list
        每个依赖项为一行，以列表的形式返回
    """

    return get_experiment_requirements(experiment_id)


# 修改实验可见性
@router.patch("/{experiment_id}/show")
async def _(experiment_id: int, request: Request):
    """修改实验是否可见

    Parameters
    ----------
    experiment_id : int
        实验id
    request : Request
        请求体，包含 show 字段

    Returns
    -------
    experiment : dict
        当前实验信息
    """

    show = (await request.json()).get("show")
    if show is None:
        return PARAMS_ERROR_422("Request parameter 'show'")
    return change_experiment_visibility(experiment_id, show)
