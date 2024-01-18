#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-17 14:56:38
@File: swanlab\db\modules\experiment.py
@IDE: vscode
@Description:
    实验数据表
"""
from ..settings import swandb
from ..model import SwanModel
from peewee import ForeignKeyField, CharField, IntegerField, TextField
from .projects import Project, DEFAULT_PROJECT_ID
from ...utils.time import create_time


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

    class Meta:
        database = swandb
        # 通过meta规定name和project_id的唯一性
        indexes = ((("name", "project_id"), True), (("index", "project_id"), True))

    id = IntegerField(primary_key=True)
    """实验id"""
    project_id = ForeignKeyField(Project, backref="experiments", default=1)
    """外键，项目id，可通过此外键反向查询项目下的所有实验"""

    run_id = IntegerField(unique=True)
    """运行时id，用于区分不同的实验保存的文件夹名称"""

    name = CharField(max_length=100, null=False)
    """实验名称，通过meta规定name和project_id的唯一性"""

    description = CharField(max_length=255)
    """实验描述"""

    index = IntegerField()
    """实验索引，用于排序，通过meta规定index和project_id的唯一性"""

    status = IntegerField(choices=[-1, 0, 1], default=0)
    """实验状态，-1: crushed, 0: running, 1: finished，在创建一个实验时，状态默认为0"""

    show = IntegerField(default=1, choices=[0, 1])
    """实验可见性，0: 不可见，1: 可见"""

    more = TextField(default="")
    """更多信息配置，json格式，将在表函数中检查并解析"""

    create_time = CharField(max_length=30, null=False)
    """创建时间"""
    update_time = CharField(max_length=30, null=False)
    """更新的时间"""

    @classmethod
    def create(
        cls,
        name: str,
        run_id: str,
        description: str = "",
        project_id: int = 1,
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
        more : dict, optional
            更多配置，在函数内部会转换为字符串格式

        Returns
        -------
        Experiment
            新建的实验实例
        """
        current_time = create_time()
        # 调用父类的create方法创建实验实例
        experiment = super().create(
            name=name,
            run_id=run_id,
            project_id=project_id,
            description=description,
            more=cls.dict2json(more),
            create_time=current_time,
            update_time=current_time,
        )

        return experiment

    @classmethod
    def create_experiment(
        cls,
        run_id,
        project_id=1,
        name="Experiment",
        description="This is a sample experiment.",
        index=0,
        status=0,
        show=1,
        more="",
    ):
        """创建实验表
        其中需要注意：run_id 必须传递！
        TODO: 是否需要在创建实验数据之前 init Project 表？
        答：不需要
        """

        # 确保 run_id 唯一
        if not cls.check_run_id(run_id):
            raise ValueError("run_id already exists")
        # 获取之前的实验数量
        experiment_count = cls.select().count()
        experiment = cls.create(
            project_id=project_id,
            run_id=run_id,
            name=name + f"-{experiment_count}",
            description=description,
            index=index,
            status=status,
            show=show,
            more=more,
            create_time=create_time(),
            update_time=create_time(),
        )

        if experiment:
            Project.update_sum(type="increase")
            return experiment
        else:
            raise Exception("Failed to create experiment")

    @classmethod
    def check_run_id(cls, run_id):
        """检查run_id是否可用"""

        return cls.filter(cls.run_id == run_id).count() == 0

    @classmethod
    def get_experiment(cls, id):
        """根据id获取实验"""

        result = cls.filter(cls.id == id)
        return SwanModel.search2dict(result)

    @classmethod
    def get_experiment_by_runid(cls, run_id):
        """根据run_id获取实验"""

        return SwanModel.search2dict(cls.filter(cls.run_id == run_id))

    @classmethod
    def get_experiment_by_name(cls, name):
        """根据name获取实验"""

        return SwanModel.search2dict(cls.filter(cls.name == name))

    @classmethod
    def get_experiments(cls):
        """获取所有的实验"""

        return SwanModel.search2dict(cls.select())

    @classmethod
    def get_tags(cls, id: int):
        """获取实验下所有的标签数据

        Parameters
        ----------
        id : int
            实验的 id

        Returns
        -------
        list
            标签列表
        """

        return SwanModel.search2list(cls.filter(cls.id == id).first().tags)

    @classmethod
    def get_charts(cls, id: int):
        """获取实验下所有的图标数据

        Parameters
        ----------
        id : int
            实验的 id

        Returns
        -------
        list
            图表列表
        """

        return SwanModel.search2list(cls.filter(cls.id == id).first().charts)

    @classmethod
    def delete_experiment(cls, id):
        """删除实验"""

        experiment = cls.filter(cls.id == id).first()
        if experiment:
            experiment.delete_instance(recursive=True)
            Project.update_sum(type="decrease")
            return True
        else:
            return False

    @classmethod
    def update_info(cls, id: int, name: str, description: str = ""):
        """更新实验信息

        Parameters
        ----------
        id : int
        name : str
        description : str, optional
        """

        experiment = cls.filter(cls.id == id).first()
        if experiment:
            experiment.name = name
            experiment.description = description
            experiment.save()
            return True
        else:
            raise Exception("Failed to update experiment")

    @classmethod
    def update_updatetime(cls, id: int):
        """实验有更新，刷新 update_time

        Parameters
        ----------
        id : int
        """

        experiment = cls.filter(cls.id == id).first()
        if experiment:
            experiment.update_time = create_time()
            experiment.save()
            return True
        else:
            raise Exception("Failed to renew update time")

    @classmethod
    def update_status(cls, id: int, status: int):
        """更新实验状态

        Parameters
        ----------
        id : int
            实验 id
        status : int
            * -1: crushed
            *  0: running
            *  1: finished
        """

        experiment = cls.filter(cls.id == id).first()
        if experiment:
            experiment.status = status
            experiment.save()
            return True
        else:
            raise Exception("Failed to update experiment status")
