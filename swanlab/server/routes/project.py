from fastapi import APIRouter

from ..controller.project import (
    get_experiments_list,
    get_project_summary,
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
    return get_project_summary()
