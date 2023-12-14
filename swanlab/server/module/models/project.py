from ....utils import lock_file
from io import TextIOWrapper
from ....env import SWANLAB_LOGS_FOLDER
import os
import ujson


class PT(object):
    """后端层面上的项目管理类，适配后端的项目管理接口，提供项目管理的相关功能"""

    path = os.path.join(SWANLAB_LOGS_FOLDER, "project.json")

    @lock_file(file_path=path, mode="r")
    def get(self, file: TextIOWrapper):
        """获取实验信息"""
        return ujson.load(file)
