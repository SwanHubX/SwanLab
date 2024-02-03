#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-02-03 16:00:13
@File: swanlab\db\utils\chart.py
@IDE: vscode
@Description:
    关于表格的一些特殊处理
    这里的函数可能在 server、data 中都有作用
"""

from peewee import fn
from ...utils import create_time
from ..models import (
    Project,
    Namespace,
    Tag,
    Source,
    Chart,
    Display,
)

from ...db import connect, NotExistedError

__to_list = Project.search2list


def transform_to_multi_exp_charts(project_id: int):
    """兼容以前没有多实验对比数据的情况
    根据早期版本的数据创建多实验对比图表

    Parameters
    ----------
    project_id : int
        project 表 id

    Raises
    ------
    NotExistedError
        project 不存在
    NotExistedError
        没有满足多实验图表的 tag
    """

    project = Project.filter(Project.id == project_id)
    if project.count() == 0:
        raise NotExistedError("Target project does not exist")

    if project.first().charts == 1:
        return "have multiple charts"

    # ---------------------------------- 未生成多实验图表 ----------------------------------

    # 1. 查找
    # - 找到所有满足条件的 tag : tag 名在多个实验中出现
    # - tag 名一样，但类型不一样，以第一个被找到的 tag 类型为基准
    # 2. 生成
    # - 生成 namespace
    # - 生成 chart
    # - 添加 chart 和 namespace 到 display
    # - 生成 source
    # 3. 标记
    # - 设置 project 中字段 charts 为 1，表示已经生成多实验表格

    # 获取名称在两个及以上实验中出现过的 tag，同时 type 一致
    tags_with_same_name = (
        Tag.select(Tag.name, Tag.type, fn.COUNT(Tag.id).alias("count"))
        .group_by(Tag.name, Tag.type)
        .having(fn.COUNT(Tag.id) > 1)
    )
    tags = [{"name": tag.name, "type": tag.type} for tag in tags_with_same_name]
    # 如果没有满足多实验图表生成条件的 tag，返回404
    if len(tags) == 0:
        raise NotExistedError("No mutiple experment charts found")
    result_lists = {}
    # 遍历查询，获取满足多实验图表生成条件的 tag
    for tag in tags:
        # 如果当前 tag 名已经存在收集列表，说明有同名不同类的 tag 都满足多图表条件，目前只取第一次出现的 type，即后续有同名不同类的 tag 都忽略
        if tag["name"] in result_lists:
            continue
        # 获取当前name的所有数据行
        target = Tag.filter(Tag.name == tag["name"], Tag.type == tag["type"])
        # 将数据行添加到结果列表
        result_lists[tag["name"]] = [
            {"experiment_name": item["experiment_id"]["name"], "tag_id": item["id"], "type": tag["type"]}
            for item in __to_list(target)
        ]
    # 至此 result_list 数据结构为：{ tag_name: [{experiment_name, tag_id, tag_type}] }，找到了所有可以生成多实验图表的 tag
    # 查找过滤后，需要生成多实验图表，生成操作需要使用原子操作
    db = connect()
    with db.atomic():
        # 1. 生成 namespace
        namespace = Namespace.create("default", project_id=project_id)
        for key in result_lists:
            # 2. 生成图表
            chart = Chart.create(key, experiment_id=None, project_id=project_id, system=0)
            # 3. 添加 chart 和 namespace 到 display
            Display.create(chart_id=chart.id, namespace_id=namespace.id)
            # 4. 生成source
            time = create_time()
            sources = [
                {
                    "tag_id": item["tag_id"],
                    "chart_id": chart.id,
                    "sort": index,
                    "create_time": time,
                    "update_time": time,
                }
                for index, item in enumerate(result_lists[key])
            ]
            Source.insert_many(sources).execute()
        # 最后将 project 表中的 charts 标注为 1，即已经生成多实验表格
        Project.update(charts=1).where(Project.id == project_id).execute()
    db.commit()
