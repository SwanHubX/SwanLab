"""
@author: cunyue
@file: mounter.py
@time: 2025/7/17 13:46
@description: 挂载项目和实验，并在 run_store 中记录相关信息
一般情况下，他与 Porter 一起使用，作为 Porter 的前置操作
"""

from typing import Optional

from swanlab.core_python import get_client, Client
from swanlab.data import namer as N
from swanlab.data.run.config import SwanLabConfig
from swanlab.data.run.metadata.hardware import is_system_key
from swanlab.data.store import RunStore, get_run_store, RemoteMetric


class Mounter:
    """
    Mounter is a context manager that mounts a project and an experiment to the run_store.
    It initializes the run_store with the project and experiment information, and provides methods to execute the mounting process.
    """

    def __init__(self):
        """
        Initializes the Mounter.
        example usage:
        >>> with Mounter() as mounter:
        >>>     mounter.execute()  # create project and experiment, and mount them to run_store
        """
        self._run_store: Optional[RunStore] = None
        self._client: Optional[Client] = None

    @property
    def run_store(self) -> RunStore:
        assert self._run_store is not None, "You can only access run_store inside the context manager."
        return self._run_store

    def __enter__(self):
        try:
            self._client = get_client()
        except ValueError:
            raise ValueError("You must create a client before using Mounter.")
        self._run_store = get_run_store()
        return self

    def execute(self):
        """
        执行挂载
        """
        run_store = self.run_store
        http = self._client
        http.mount_project(run_store.project, run_store.workspace, run_store.visibility)
        # 1. 挂载实验前执行必要的参数生成操作
        exp_count = http.history_exp_count
        run_store.run_name = N.generate_name(exp_count) if run_store.run_name is None else run_store.run_name
        run_store.description = "" if run_store.description is None else run_store.description
        run_store.run_colors = N.generate_colors(exp_count) if run_store.run_colors is None else run_store.run_colors
        run_store.tags = [] if run_store.tags is None else run_store.tags
        # 2. 挂载实验
        try:
            new = http.mount_exp(
                exp_name=run_store.run_name,
                colors=run_store.run_colors,
                description=run_store.description,
                tags=run_store.tags,
                cuid=run_store.run_id,
                must_exist=run_store.resume == 'must',
            )
            run_store.new = new
        except RuntimeError:
            raise RuntimeError(
                "When resume=must, the experiment must exist in project {}. Please check your parameters.".format(
                    http.proj.name
                )
            )
        run_store.run_id = http.exp_id
        # 如果不是新实验，需要获取最新的实验配置和指标数据进行本地同步
        if not new:
            # 1. 解析 config，保存到 run_store.config
            config = http.exp.config
            run_store.config = SwanLabConfig.revert_config(config or {})
            # 2. 获取最新的指标数据，解析为 run_store.metrics
            # 指标总结数据，{log:[{key: '', step: int}] or None, media: [{key: str, step: int}, ...] or None, scalar: [{key: str, step: int}, ...] or None}
            summaries, _ = http.get(
                "/house/metrics/summaries/{}/{}".format(
                    http.exp.root_proj_cuid or http.proj.cuid, http.exp.root_exp_cuid or http.exp_id
                ),
                {"all": True},
            )
            # 设置历史指标条数
            run_store.log_epoch = summaries.get("log", [{}])[0].get("step", 0)
            # 列响应
            columns_resp, _ = http.get("/experiment/{}/column".format(http.exp_id), {"all": True})
            # 列信息, [{error: dict or None, key:str, type:str, class:str}]
            columns = columns_resp.get("list", [])
            # key -> (column_type, column_class, error, latest step)
            metrics: RemoteMetric = {}
            for column in columns:
                # 从列信息中获取指标信息
                key = column["key"]
                column_type = column["type"]
                column_class = column["class"]
                error = column.get("error", None)
                if column_class == "SYSTEM" and not is_system_key(key):
                    # 只记录 sdk 生成的系统指标
                    continue
                # 从总结数据中获取最新的 step
                # 这里需要同时查找 media 和 scalar
                latest_step = None
                for scalar_summary in summaries.get("scalar") or []:
                    if scalar_summary["key"] == key:
                        latest_step = scalar_summary["step"]
                        break
                if latest_step is None:
                    for media_summary in summaries.get("media") or []:
                        if media_summary["key"] == key:
                            latest_step = media_summary["step"]
                            break
                metrics[key] = (column_type, column_class, error, latest_step)
            run_store.metrics = metrics

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._run_store = None
