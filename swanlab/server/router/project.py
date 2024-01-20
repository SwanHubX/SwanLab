from fastapi import APIRouter, Request

from ..controller.project import (
    # 列出实验列表
    get_experiments_list,
    # 获取项目总结信息
    get_project_summary,
    # 修改项目信息
    update_project_info,
    # 删除项目
    delete_project,
)

router = APIRouter()


@router.get("")
def _():
    """获取项目下的实验列表

    Returns
    -------
    list[dict]:
        项目下的实验列表
    """

    return get_experiments_list()


@router.get("/summaries")
def _():
    """获取项目总结信息

    Returns
    -------
    dict
        项目下每个实验的总结信息
    """

    return get_project_summary()


@router.patch("/update")
async def _(request: Request, project_id: int = None):
    """修改项目信息

    Parameters
    ----------
    request : Request
        请求体，应该包含 name 和 description 字段
    project_id : int, optional
        by default None

    Returns
    -------
    dict
        修改后完整的 Project 信息
    """

    if project_id is None:
        return await update_project_info(request)
    return await update_project_info(request, project_id)


@router.delete("")
async def _(project_id: int = None):
    """删除项目

    Parameters
    ----------
    project_id : int, optional
        默认为第一个项目
    """

    if project_id is None:
        return await delete_project()
    return await delete_project(project_id)
