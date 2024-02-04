from fastapi import APIRouter, Request

from ..controller.project import (
    # 列出实验列表
    get_project_info,
    # 获取项目总结信息
    get_project_summary,
    # 修改项目信息
    update_project_info,
    # 删除项目
    delete_project,
    # 获取多实验对比图表
    get_project_charts,
    # 默认项目id
    DEFAULT_PROJECT_ID,
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

    return get_project_info()


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
        修改后的项目信息
        - name
        - description
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


@router.get("/charts")
async def _(project_id: int = DEFAULT_PROJECT_ID):
    """获取多实验对比图表数据,并且考虑往期版本兼容性
    1. 如果当前项目的chart字段为0，先生成多实验对比数据，跳转步骤2
    2. 依据规则获取所有实验的图表数据
    """

    return await get_project_charts(project_id)
