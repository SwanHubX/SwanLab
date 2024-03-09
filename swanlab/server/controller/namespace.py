import os
import ujson
import shutil
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409, NOT_FOUND_404
from fastapi import Request

import yaml
from typing import List, Dict

from ...db import connect, NotExistedError

from ...db import (
    Namespace,
)


__to_list = Namespace.search2list


# 修改namespace可见性
def change_namespace_visibility(namespace_id: int, opened: int):
    """修改namespace可见性

    Parameters
    ----------
    namespace_id : int
        _description_
    opened : int
        在一个实验下，该namespace是否可见，true -> 1 , false -> 0

    Returns
    -------
    _type_:dict
        当前namespace信息
    """
    try:
        namespace: Namespace = Namespace.get_by_id(namespace_id)
    except NotExistedError:
        return NOT_FOUND_404("Namespace with id {} does not exist.".format(namespace_id))
    namespace.opened = 1 if opened else 0
    namespace.save()
    return SUCCESS_200(None)
