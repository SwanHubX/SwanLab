from fastapi import APIRouter
from ..module.resp import PARAMS_ERROR_422

from ..controller.namespace import (
    # 修改namespace可见性
    change_namespace_visibility,
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

    opened = (await request.json()).get("opened")
    if opened is None:
        return PARAMS_ERROR_422("Request parameter 'show'")
    return change_namespace_visibility(namespace_id, opened)
