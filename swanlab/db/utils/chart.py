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
            chart = Chart.create(
                key, experiment_id=None, project_id=project_id, system=0, type=result_lists[key][0]["type"]
            )
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


def add_multi_chart(project_id: int, experiment_id: int, tag_id: int, chart_id: int):
    """添加新的多实验对比图表
    1. 当前 tag 已有多实验图表，将该 tag 添加 source 即可
    2. 当前 tag 满足条件但未创建多实验图表，进行创建和添加

    Parameters
    ----------
    project_id : int
        项目 id
    experiment_id : int
        实验 id
    tag_id : int
        tag id
    chart_id : int
        tag 伴生 chart 的 id
    """

    try:
        project = Project.filter(Project.id == project_id).first()
    except Exception:
        raise NotExistedError("Target project does not exist")

    # ---------------------------------- charts 为 0 说明还没有生成过多实验图表 ----------------------------------
    if project.charts == 0:
        try:
            transform_to_multi_exp_charts(project_id)
        except NotExistedError:
            # 转化失败（没有可转换的）
            return False
        return True

    # ---------------------------------- charts 为 1 说明已经生成过多实验图表 ----------------------------------
    tag = Tag.get_by_id(tag_id)
    # 找出所有和该 tag 同名同类的行
    tags = Tag.filter(Tag.name == tag.name, Tag.type == tag.type)

    # 检查该 tag 是否已经有满足的图表了 => 下面的逻辑只考虑了单 namesapce

    if tags.count() <= 1:
        # 如果只有一个，说明当前 tag 是第一个，不满足添加和创建条件，略过
        return "No tag to create multiple experimental chart"

    # 起码三个
    if tags.count() > 2:
        # 当前 tag 是否有对应的图表，因为可能同名不同类，当不同类都满足条件时，后出现的类型忽略，只有最初的类型才可绑定到对应表
        charts = Chart.filter(Chart.project_id == project_id, Chart.name == tag.name, Chart.type == tag.type)
        if charts.count() == 0:
            # 之前有同名不同类的 tag 已经创建了多实验对比图表
            return "Previously, there were tags with the same name but different types, and the chart has been created"
        # 说明之前已经有这种表格了，且当前 tag 不是第一个，现在只需将 tag 插入多实验图表
        chart = charts.first()
        # 将 tag 和相应的多实验对比图表对应起来 => 使用 source
        Source.create(tag.id, chart.id)
        return True

    # 如果有两个，说明该 tag 在添加后，刚好满足图表创建条件
    # 但是这个时候需要注意，是否有同名不同类的 tag 已经设置了多实验图表
    charts = Chart.filter(Chart.project_id == project_id, Chart.name == tag.name)
    if not charts.count() == 0:
        c = charts.first()
        return f"There are different types of tags called this name, and the chart named {c.name} has inited in type '{c.type}'"

    # 此时刚好有两个同名同类的 tag，执行创建新多实验对比图表操作：
    # - 检查(生成) namespace
    # - 生成 chart
    # - 添加 chart 和 namespace 到 display
    # - 生成 source

    db = connect()
    with db.atomic():
        # 1. 获取 namespace，此时 project.charts 为 1，project 的默认 namsepace 一定存在
        namespace = Namespace.filter(Namespace.project_id == project_id).first()
        # 2. 生成图表
        accompanying_chart = Chart.get_by_id(chart_id)
        # 通过 tag 的伴生图表获取部分信息
        chart = Chart.create(
            name=accompanying_chart.name,
            description=accompanying_chart.description,
            project_id=project_id,
            system=0,
            type=accompanying_chart.type,
            reference=accompanying_chart.reference,
            config=accompanying_chart.config,
            more=accompanying_chart.more,
        )
        # 3. 添加 chart 和 namespace 到 display
        Display.create(chart.id, namespace.id)
        # 4. 生成 source
        time = create_time()
        sources = [
            {
                "tag_id": tag.id,
                "chart_id": chart.id,
                "sort": index,
                "create_time": time,
                "update_time": time,
            }
            for index, tag in enumerate(tags)
        ]
        Source.insert_many(sources).execute()
    db.commit()
    return True
