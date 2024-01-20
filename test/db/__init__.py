#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-20 14:23:18
@File: test\db\__init__.py
@IDE: vscode
@Description:
    每次一个个添加数据太麻烦了，写一个统筹的东西
"""

from swanlab.db import connect
from swanlab.db import (
    Project,
    Experiment,
    Chart,
    Tag,
    Namespace,
    Display,
    Source,
)

swandb = connect()

with swandb.atomic():
    Project.init()

    Experiment.create("experiment1", "run1")
    Experiment.create("experiment2", "run2")
    Experiment.create("experiment3", "run3")

    Tag.create(experiment_id=1, name="tag1", type="int")
    Tag.create(experiment_id=2, name="tag2", type="int")
    Tag.create(experiment_id=3, name="tag3", type="int")
    Tag.create(experiment_id=1, name="tag11", type="int")
    Tag.create(experiment_id=2, name="tag1", type="int")

    Chart.create(name="chart1", description="123", project_id=1)
    Chart.create(name="chart2", description="123", project_id=1)
    Chart.create(name="chart11", description="123", experiment_id=1)
    Chart.create(name="chart22", description="123", experiment_id=1)

    Namespace.create(name="name1", experiment_id=1)
    Namespace.create(name="name2", experiment_id=1)
    Namespace.create(name="name3", experiment_id=1)
    Namespace.create(name="name11", project_id=1)
    Namespace.create(name="name22", project_id=1)
