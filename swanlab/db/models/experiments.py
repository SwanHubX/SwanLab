#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 14:56:38
@File: swanlab\db\modules\experiment.py
@IDE: vscode
@Description:
    实验数据表
"""
import os.path

from ..model import SwanModel
from peewee import ForeignKeyField, CharField, IntegerField, TextField, IntegrityError, Check, DatabaseProxy
from .projects import Project
from ..error import ExistedError, NotExistedError, ForeignProNotExistedError
from swanlab.package import get_package_version
from ...utils import generate_color
from ...utils import create_time
from swanlab.env import get_swanlog_dir
import shutil


# 定义模型类
class Experiment(SwanModel):
    """实验表
    """

    # 实验运行时状态符
    RUNNING_STATUS = 0
    # 实验停止时状态符
    STOPPED_STATUS = -1
    # 实验结束时状态符
    FINISHED_STATUS = 1

    class Meta:
        database = DatabaseProxy()
        # 通过meta规定name和project_id的唯一性
        indexes = ((("name", "project_id"), True), (("sort", "project_id"), True))
        # sort必须大于等于0
        constraints = [Check("sort >= 0")]

    id = IntegerField(primary_key=True)
    """实验id"""
    project_id = ForeignKeyField(Project, backref="experiments", on_delete="CASCADE", on_update="CASCADE", default=1)
    """外键，项目id，可通过此外键反向查询项目下的所有实验"""

    run_id = TextField(unique=True)
    """运行时id，用于区分不同的实验保存的文件夹名称"""

    name = CharField(max_length=100, null=False)
    """实验名称，通过meta规定name和project_id的唯一性"""

    description = CharField(max_length=255, null=True)
    """实验描述"""

    sort = IntegerField(null=False)
    """实验索引，用于排序，通过meta规定index和project_id的唯一性"""

    status = IntegerField(choices=[-1, 0, 1], default=0)
    """实验状态，-1: crushed, 0: running, 1: finished，在创建一个实验时，状态默认为0"""

    show = IntegerField(default=1, choices=[0, 1])
    """实验可见性，0: 不可见，1: 可见"""

    light = CharField(max_length=20, null=True)
    """亮色主题颜色"""
    dark = CharField(max_length=20, null=True)
    """暗色主题颜色"""

    pinned_opened = IntegerField(default=1, choices=[0, 1])
    """实验图表置顶部分是否打开，默认值为1，表示打开"""

    hidden_opened = IntegerField(default=0, choices=[0, 1])
    """实验图表隐藏部分是否打开，默认值为0,表示关闭"""

    more = TextField(null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""

    version = CharField(max_length=30, null=False)
    """创建时的版本号"""

    create_time = CharField(max_length=30, null=False)
    """创建时间"""

    finish_time = CharField(max_length=30, null=True, default=None)
    """实验结束时间"""

    update_time = CharField(max_length=30, null=False)
    """更新的时间"""

    def __dict__(self):
        return {
            "id": self.id,
            "project_id": self.project_id.__dict__(),
            "run_id": self.run_id,
            "name": self.name,
            "description": self.description,
            "sort": self.sort,
            "status": self.status,
            "show": self.show,
            "more": self.more,
            "version": self.version,
            "create_time": self.create_time,
            "update_time": self.update_time,
            "finish_time": self.finish_time,
        }

    @classmethod
    def create(
        cls,
        name: str,
        run_id: str,
        description: str = None,
        project_id: int = Project.DEFAULT_PROJECT_ID,
        num: int = None,
        more: dict = None,
    ) -> "Experiment":
        """覆写继承的create方法，创建一个新的实验

        Parameters
        ----------
        name : str
            实验名称
        run_id : str
            运行时信息
        description : str
            实验描述，默认为空字符串
        project_id : int, optional
            关联的项目id，由于目前为单项目模式，所以不需要手动设置此字段，默认置为DEFAULT_PROJECT_ID
        num : int, optional
            历史实验数量，如果不指定则自动计算
        more : dict, optional
            更多配置，在函数内部会转换为字符串格式

        Returns
        -------
        Experiment
            新建的实验实例

        Raises
        -------
        ForeignProNotExistedError
            项目不存在
        ExistedError
            实验已经存在
        """
        # 检查项目是否存在
        if not Project.select().where(Project.id == project_id).exists():
            raise ForeignProNotExistedError("项目不存在")
        # 这个sum是+1以后的值
        _sum = Project.increase_sum(project_id)
        # 调用父类的create方法创建实验实例
        light, dark = generate_color(_sum if num is None else num + 1)
        try:
            return super().create(
                name=name,
                run_id=run_id,
                project_id=project_id,
                description=description,
                more=cls.dict2json(more),
                sort=_sum,
                version=get_package_version(),
                light=light,
                dark=dark,
            )
        except IntegrityError:
            raise ExistedError("实验已经存在")

    @classmethod
    def get(cls, *args, **kwargs) -> "Experiment":
        try:
            return super().get(*args, **kwargs)
        except:
            raise NotExistedError("Experiment does not exist: {}".format(kwargs))

    def update_status(self, status: int):
        """更新实验状态

        Parameters
        ----------
        status : int
            实验状态
        """
        self.status = status
        self.finish_time = create_time()
        self.save()

    @classmethod
    def purely_delete(cls, run_id: str):
        """
        删除某一个实验，顺便删除实验的文件夹

        Parameters
        ----------
        run_id : int
            实验run_id
        """
        # 获取实验实例
        try:
            exp = cls.get(run_id=run_id)
            # 删除实验
            exp.delete_instance()
            # 更新项目的实验数量
            Project.decrease_sum()
            # 删除实验文件夹
            shutil.rmtree(os.path.join(get_swanlog_dir(), exp.run_id.__str__()))
        except NotExistedError:
            shutil.rmtree(os.path.join(get_swanlog_dir(), run_id))
