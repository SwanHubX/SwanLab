"""
@author: cunyue
@file: sync_utils.py
@time: 2025/7/21 17:20
@description: swanlab sync 函数的工具函数
"""

from typing import Optional, Union, Literal

from ..data.store import RunStore
from ..proto.v0 import Project, Experiment


def set_run_store(
    run_store: RunStore,
    proj: Project,
    exp: Experiment,
    project: Optional[str] = None,
    workspace: Optional[str] = None,
    id: Union[str, Literal['auto', 'new']] = "new",
):
    """
    在同步之前设置运行存储的相关信息。
    :param run_store: RunStore 实例
    :param proj: 项目信息
    :param exp: 实验信息
    :param project: 指定的项目名称，如果为 None 则使用 proj.name
    :param workspace: 指定的工作空间名称，如果为 None 则使用 proj.workspace
    :param id: 实验 ID，可以是 'new', 'auto' 或具体的 ID 字符串 本函数不做字符串格式检查
    """
    # 设置项目信息
    run_store.project = project or proj.name
    run_store.workspace = workspace or proj.workspace
    run_store.visibility = proj.public
    # 设置实验信息
    run_store.run_name = exp.name
    run_store.description = exp.description
    run_store.tags = exp.tags
    run_store.run_colors = (exp.colors[0], exp.colors[1])
    # 设置实验 id 和 resume 模式
    # a. id 为 new，则 resume 为 never, id 为 None
    # b. id 为 auto，则 resume 为 allow, id 为 exp.id
    # c. 其他情况，则 resume 为 must, id 为 id
    if id == "new":
        run_store.resume = 'never'
        run_store.run_id = None
    elif id == "auto":
        run_store.resume = 'allow'
        run_store.run_id = exp.id
    else:
        run_store.resume = 'must'
        run_store.run_id = id
