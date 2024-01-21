#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 14:56:38
@File: swanlab\db\modules\experiment.py
@IDE: vscode
@Description:
    实验数据表
"""
from ..model import SwanModel
from peewee import ForeignKeyField, CharField, IntegerField, TextField, IntegrityError, Check, DatabaseProxy
from .projects import Project
from ..error import ExistedError, NotExistedError
from ...utils.time import create_time
from ...utils.package import get_package_version
from ...utils import generate_color


# 定义模型类
class Experiment(SwanModel):
    """实验表

    Attributes
    ----------
    tags: list of Tag
        由 Tag 表中外键反链接生成的tag数据列表
    charts: list of Chart
        由 Chart 表中外键反链接生成的chart数据列表
    namespaces: list of Namespace
        由 Namespace 表中外键反链接生成的命名空间数据列表
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

    more = TextField(null=True)
    """更多信息配置，json格式，将在表函数中检查并解析"""

    version = CharField(max_length=30, null=False)
    """创建时的版本号"""

    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新的时间"""

    def __dict__(self):
        return {
            "id": self.id,
            "project_id": self.project_id,
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
        }

    @classmethod
    def create(
        cls,
        name: str,
        run_id: str,
        description: str = None,
        project_id: int = Project.DEFAULT_PROJECT_ID,
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
        sort : int
            实验索引，默认为None，如果为None，则会自动设置为当前项目下的实验数量
        project_id : int, optional
            关联的项目id，由于目前为单项目模式，所以不需要手动设置此字段，默认置为DEFAULT_PROJECT_ID
        more : dict, optional
            更多配置，在函数内部会转换为字符串格式

        Returns
        -------
        Experiment
            新建的实验实例

        Raises
        -------
        NotExistedError
            项目不存在
        ExistedError
            实验已经存在
        """
        current_time = create_time()
        # 检查项目是否存在
        if not Project.select().where(Project.id == project_id).exists():
            raise NotExistedError("项目不存在")

        # 这个sum是+1以后的值，所以需要-1
        sum = Project.increase_sum(project_id)
        # 自动设置索引，为当前项目下的实验数量
        sort = cls.select().where(cls.project_id == project_id).count()
        # 调用父类的create方法创建实验实例
        try:
            return super().create(
                name=name,
                run_id=run_id,
                project_id=project_id,
                description=description,
                more=cls.dict2json(more),
                sort=sort,
                version=get_package_version(),
                light=generate_color(sum),
                create_time=current_time,
                update_time=current_time,
            )
        except IntegrityError:
            raise ExistedError("实验已经存在")

    @classmethod
    def get(cls, *args, **kwargs):
        """覆写继承的get方法，通过id获取实验实例

        Parameters
        ----------
        id : int
            实验id

        Returns
        -------
        Experiment
            实验实例

        Raises
        -------
        NotExistedError
            实验不存在
        """
        try:
            return super().get(*args, **kwargs)
        except:
            raise NotExistedError("Experiment does not exist: {}".format(kwargs))
