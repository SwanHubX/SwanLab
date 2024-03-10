from fastapi import APIRouter
from ..module.resp import PARAMS_ERROR_422

from ..controller.chart import (
    update_charts_status,
)
from fastapi import Request
from urllib.parse import quote

router = APIRouter()


@router.patch("/{chart_id}/status")
async def _(chart_id: int, request: Request):
    """修改图表状态

    Parameters
    ----------
    chart_id : int
        图表id
    request : Request
        请求体，包含 status 字段
    Returns
    -------
    groups : list
        图表列表
    """

    data = await request.json()
    status = data.get("status")

    if status is None:
        return PARAMS_ERROR_422("Request parameter 'status'")
    return update_charts_status(chart_id, status)
