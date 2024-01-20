from fastapi import APIRouter

from ..controller.experiment import (
    get_experiment_info,
)

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
