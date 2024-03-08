import os
import ujson
import shutil
from ..module.resp import SUCCESS_200, DATA_ERROR_500, CONFLICT_409, NOT_FOUND_404
from fastapi import Request

import yaml
from ...log import swanlog
from typing import List, Dict

from ...db import connect, NotExistedError

from ...db import (
    Namespace,
)


__to_list = Namespace.search2list


# 修改namespace下可见性
def change_namespace_visibility(namespace_id: int, opened: bool):
    namespace = Namespace.get_by_id(namespace_id)
    namespace.opened = opened
    namespace.save()
    return SUCCESS_200({"namespace": namespace_id})
