"""
@author: cunyue
@file: utils.py
@time: 2026/5/19 15:55
@description: Core 服务共享工具函数
"""

from dataclasses import dataclass

from swanlab.proto.swanlab.run.v1.run_pb2 import StartRecord
from swanlab.sdk.internal.core_python.api.experiment import (
    create_or_resume_experiment as api_create_or_resume_experiment,
)
from swanlab.sdk.internal.core_python.api.project import get_or_create_project, get_project
from swanlab.sdk.internal.pkg import adapter
from swanlab.sdk.typings.core_python.api.experiment import InitExperimentType
from swanlab.sdk.typings.core_python.api.project import ProjectType
from swanlab.utils.experiment import generate_color, generate_name


@dataclass(frozen=True)
class PrepareExperimentStartResult:
    username: str
    project: str
    project_info: ProjectType
    experiment: InitExperimentType
    new_experiment: bool
    name: str
    color: str


def prepare_experiment_start(record: StartRecord) -> PrepareExperimentStartResult:
    """
    创建或恢复运行对应的项目和实验
    """
    project_data = get_or_create_project(
        username=record.workspace,
        name=record.project,
        public=record.public,
    )
    username, project = project_data["username"], project_data["name"]

    project_info = get_project(username=username, name=project)
    history_experiment_count = project_info["_count"]["experiments"]
    name = record.name or generate_name(history_experiment_count)
    color = record.color or generate_color(history_experiment_count)
    resume = adapter.resume[record.resume]

    experiment, new_experiment = api_create_or_resume_experiment(
        username,
        project,
        name=name,
        resume=resume,
        run_id=record.id,
        color=color,
        description=record.description,
        job_type=record.job_type,
        group=record.group,
        tags=list(record.tags),
        created_at=record.started_at,
    )
    assert experiment.get("name"), "create_or_resume_experiment() returned an experiment without a name."

    return PrepareExperimentStartResult(
        username=username,
        project=project,
        project_info=project_info,
        experiment=experiment,
        new_experiment=new_experiment,
        name=name,
        color=color,
    )
