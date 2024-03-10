from fastapi import APIRouter
from ..module.resp import PARAMS_ERROR_422

from ..controller.namespace import (
    # 修改namespace开启关闭状态
    change_namespace_opened,
)
from fastapi import Request
from urllib.parse import quote

router = APIRouter()


# 修改namespace可见性
@router.patch("/{namespace_id}/opened")
async def _(namespace_id: int, request: Request):
    """修改namespace可见性

    Parameters
    ----------
    experiment_id : int
        实验id
    request : Request
        请求体，包含 opened 字段
    Returns
    -------
    experiment : dict
        当前实验信息
    """

    data = await request.json()
    opened = data.get("opened")
    experiment_id = data.get("experiment_id")
    project_id = data.get("project_id")

    if opened is None:
        return PARAMS_ERROR_422("Request parameter 'opened'")
    if experiment_id is None and project_id is None:
        return PARAMS_ERROR_422("Request parameter 'experiment_id' or 'project_id'")
    return change_namespace_opened(namespace_id, opened, experiment_id, project_id)
