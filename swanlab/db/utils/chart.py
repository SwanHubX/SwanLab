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
from ..models import (
    Project,
    Namespace,
    Tag,
    Source,
    Chart,
    Display,
)

from ..db_connect import connect
from ..error import (
    NotExistedError,
    ChartTypeError,
)


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
    IndexError
        project 已经生成过多实验图表
    NotExistedError
        没有满足多实验图表的 tag
    """

    project = Project.filter(Project.id == project_id)
    if project.count() == 0:
        raise NotExistedError("Target project does not exist")

    if project.first().charts == 1:
        raise IndexError("project.charts must be 0, which means no multiple charts")

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

    # 根据名称和类型分类获取 tag
    tags_with_same_name = Tag.select(Tag.name, Tag.type, fn.COUNT(Tag.id).alias("count")).group_by(Tag.name, Tag.type)
    tags = [{"name": tag.name, "type": tag.type} for tag in tags_with_same_name]
    # 如果没有满足多实验图表生成条件的 tag，抛出异常
    # 这是一个非预期异常，因为按道理不会出现这种情况,这种情况也不会写到文档中，仅仅是调试方便
    if len(tags) == 0:
        raise ValueError("No mutiple experment charts found")
    # 用于收集满足多实验图表生成条件的 tag
    result_lists = {}
    # 遍历查询，获取满足多实验图表生成条件的 tag
    for tag in tags:
        # 如果当前 tag 名已经存在收集列表，说明有同名不同类的 tag 都满足多图表条件，目前只取第一次出现的 type，即后续有同名不同类的 tag 都忽略
        if tag["name"] in result_lists and not tag["type"] == Tag.filter(Tag.name == tag["name"]).first().type:
            continue
        # 获取当前name的所有数据行
        target = Tag.filter(Tag.name == tag["name"], Tag.type == tag["type"])
        # 将数据行添加到结果列表
        result_lists[tag["name"]] = [
            {
                "experiment_id": item["experiment_id"]["id"],
                "experiment_name": item["experiment_id"]["name"],
                "tag_id": item["id"],
            }
            for item in Project.search2list(target)
        ]
    # 至此 result_list 数据结构为：{ tag_name: [{experiment_id, experiment_name, tag_id}] }，找到了所有可以生成多实验图表的 tag
    # 查找过滤后，需要生成多实验图表，生成操作需要使用原子操作
    db = connect()
    with db.atomic():
        for tag_name in result_lists:
            tag_list = result_lists[tag_name]
            # 1. 获取/生成 namespace
            # 获取 tag 对应的伴生 namespace 名
            accompanying_chart = Chart.filter(
                Chart.name == tag_name, Chart.experiment_id == tag_list[0]["experiment_id"], Chart.system != 0
            ).first()
            ns_name = Display.filter(Display.chart_id == accompanying_chart.id).first().namespace_id.name
            # 如果属于项目的同名 namespace 已经存在，直接使用这个
            namespaces = Namespace.filter(Namespace.name == ns_name, Namespace.project_id == project_id)
            if namespaces.count():
                namespace = namespaces.first()
            # 没有该命名空间，需创建
            else:
                namespace = Namespace.create(ns_name, project_id=project_id)
            # 2. 生成图表
            chart_type = Source.filter(Source.tag_id == tag_list[0]["tag_id"]).first().chart_id.type
            chart = Chart.create(
                tag_name,
                experiment_id=None,
                project_id=project_id,
                system=0,
                type=chart_type,
            )
            # 3. 添加 chart 和 namespace 到 display
            Display.create(chart_id=chart.id, namespace_id=namespace.id)
            # 4. 生成source
            sources = [
                {
                    "tag_id": item["tag_id"],
                    "chart_id": chart.id,
                    "sort": index,
                    "error": Source.filter(Source.tag_id == item["tag_id"]).first().error,
                }
                for index, item in enumerate(tag_list)
            ]
            Source.insert_many(sources).execute()
        # 最后将 project 表中的 charts 标注为 1，即已经生成多实验表格
        Project.update(charts=1).where(Project.id == project_id).execute()
    db.commit()


def add_multi_chart(
    tag_id: int,
    chart_id: int,
    project_id: int = None,
):
    """添加新的多实验对比图表
    1. 当前 tag 已有多实验图表，将该 tag 添加 source 即可
    2. 当前 tag 满足条件但未创建多实验图表，进行创建和添加

    Parameters
    ----------
    tag_id : int
        tag id
    chart_id : int
        tag 伴生 chart 的 id
    project_id : int
        项目 id, 由于目前为单项目，所以可以不传
    Raises
    ------
    ChartTypeError
        当前 tag 类型和期望 chart 类型不一致
    IndexError
        project 已经生成过多实验图表, 无法再次生成
    """
    if project_id is None:
        project_id = Project.DEFAULT_PROJECT_ID

    try:
        project = Project.filter(Project.id == project_id).first()
    except Exception:
        raise NotExistedError("Target project does not exist")

    # ---------------------------------- charts 为 0 说明还没有生成过多实验图表 ----------------------------------

    if project.charts == 0:
        return transform_to_multi_exp_charts(project_id)

    # ---------------------------------- charts 为 1 说明已经生成过多实验图表 ----------------------------------
    # 当前新增 tag
    tag = Tag.get_by_id(tag_id)
    # 属于项目的所有与 tag 同名的 chart，用于判断是否已经生成该 tag 的多实验对比图
    charts = Chart.filter(Chart.project_id == project_id, Chart.name == tag.name)
    chart = charts.first()

    # 已经有对应名称的 chart
    if charts.count() != 0:
        # 获取图所属 tag 的类型
        if chart.sources.first().tag_id.type != tag.type:
            # 同名 chart 已生成，但是类型不满足
            raise ChartTypeError("Error tag type")
        else:
            # 如果该 chart 和 tag 同类，直接添加
            return Source.create(tag_id, chart.id, error=Source.filter(Source.tag_id == tag_id).first().error)

    # 如果没有对应名称的 chart, 执行创建新多实验对比图表操作
    # 此时该 tag 应该是第一次出现
    # 获取/创建 namespace
    # - 生成或者查找 namespace
    # - 生成 chart
    # - 添加 chart 和 namespace 到 display
    # - 生成 source

    db = connect()
    with db.atomic():
        # 1. 获取 namespace
        name = Display.filter(Display.chart_id == chart_id).first().namespace_id.name
        namespaces = Namespace.filter(Namespace.project_id == project_id, Namespace.name == name)
        # 如果属于项目的同名 namespace 已经存在，直接使用这个
        if Namespace.filter(Namespace.project_id == project_id, Namespace.name == name).count() > 0:
            namespace: Namespace = namespaces.first()
        # 没有该命名空间，需创建, 生成 namespace
        else:
            namespace: Namespace = Namespace.create(name, project_id=project_id)
        # 2. 生成图表
        accompanying_chart: Chart = Chart.get_by_id(chart_id)
        # 通过 tag 的伴生图表获取部分信息
        chart = Chart.create(
            name=accompanying_chart.name,
            project_id=project_id,
            # 不是创建tag时自动生成的图表
            system=0,
            type=accompanying_chart.type,
            reference=accompanying_chart.reference,
        )
        # 3. 添加 chart 和 namespace 到 display
        Display.create(chart.id, namespace.id)
        # 4. 生成 source
        Source.create(
            tag_id,
            chart.id,
            error=Source.filter(Source.tag_id == tag_id).first().error,
        )
    db.commit()
