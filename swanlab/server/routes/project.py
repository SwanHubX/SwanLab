from fastapi import APIRouter

from ..controller.project import (
    list_experiments,
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
    return list_experiments()
