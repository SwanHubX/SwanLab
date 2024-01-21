from fastapi import APIRouter

from ..controller.experiment import (
    # 获取实验信息
    get_experiment_info,
    # 获取 tag 信息
    get_tag_data,
    # 获取实验状态
    get_experiment_status,
    # 获取实验的总结数据
    get_experiment_summary,
)

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
