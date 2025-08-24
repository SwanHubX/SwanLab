#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: experiment.py
@IDE: pycharm
@Description:
    实验相关的开放API
"""
from typing import List

from swanlab.api.base import ApiBase, ApiHTTP
from swanlab.api.types import ApiResponse, Experiment, Pagination

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = None


class ExperimentAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)

    @classmethod
    def parse(cls, body: dict) -> Experiment:
        return Experiment.model_validate({
            "cuid": body.get("cuid") or "",
            "name": body.get("name") or "",
            "description": body.get("description") or "",
            "state": body.get("state") or "",
            "show": bool(body.get("show")),
            "createdAt": body.get("createdAt") or "",
            "finishedAt": body.get("finishedAt") or "",
            "user": {
                "username": (body.get("user") or {}).get("username") or "",
                "name": (body.get("user") or {}).get("name") or "",
            },
            "profile": body.get("profile") or {}
        })

    def get_experiment(
            self,
            username: str,
            projname: str,
            exp_id: str
    ) -> ApiResponse[Experiment]:
        """
        获取实验信息

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            exp_id (str): 实验CUID
        """
        resp: ApiResponse = self.http.service.get_exp_info(username=username, project=projname, exp_id=exp_id)
        if resp.errmsg:
            return resp
        resp.data = ExperimentAPI.parse(resp.data)
        return resp

    def delete_experiment(
            self,
            username: str,
            projname: str,
            exp_id: str
    ):
        """
        删除实验

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            exp_id (str): 实验CUID
        """
        return self.http.delete(f"/project/{username}/{projname}/runs/{exp_id}", params={})

    def list_experiments(
            self,
            username: str,
            projname: str,
            page: int = 1,
            size: int = 10
    ) -> ApiResponse[Pagination[Experiment]]:
        """
        分页获取项目下的实验列表
        该接口返回的实验profile只包含config(用户自定义的配置)

        Args:
            username (str): 工作空间名
            projname (str): 项目名
            page (int): 页码, 默认为1
            size (int): 每页大小, 默认为10
        """
        resp = self.http.get(f"/project/{username}/{projname}/runs", params={"page": page, "size": size})

        if resp.errmsg:
            return resp
        exps = resp.data
        resp.data = Pagination[Experiment].model_validate({
            "total": exps.get("total", 0),
            "list": [ExperimentAPI.parse(e) for e in exps.get("list", [])]
        })
        return resp

    def get_summary(
            self,
            exp_id: str,
            pro_id: str,
            root_exp_id: str,
            root_pro_id: str
    ) -> ApiResponse[dict]:
        """
        获取实验的summary信息
        从House获取, 需要考虑克隆实验

        Args:
            exp_id (str): 实验CUID
            pro_id (str): 项目CUID
            root_exp_id (str): 根实验CUID
            root_pro_id (str): 根项目CUID
        """
        data = {
            "experimentId": exp_id,
            "projectId": pro_id,
        }
        if root_exp_id and root_pro_id:
            data["rootExpId"] = root_exp_id
            data["rootProId"] = root_pro_id

        resp = self.http.post("/house/metrics/summaries", data=[data], params={})
        if resp.errmsg:
            return resp

        resp.data = list(resp.data.values())[0]
        resp.data = {
            k: {
                "step": v.get("step"),
                "value": v.get("value"),
                "min": {
                    "step": v.get("min").get("index"),
                    "value": v.get("min").get("data"),
                },
                "max": {
                    "step": v.get("max").get("index"),
                    "value": v.get("max").get("data"),
                }
            }
            for k, v in resp.data.items()
        }
        return resp

    def get_metrics(
            self,
            exp_id: str,
            keys: List[str],
    ) -> ApiResponse[DataFrame]:
        """
        获取实验的指标数据, 可选择若干由用户自定义的列

        Args:
            exp_id (str): 实验CUID
            keys (list[str]): 指标key列表
        
        Returns:
            ApiResponse[DataFrame]:
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("OpenApi.get_metrics requires pandas module. Install with 'pip install pandas'.")

        # 去重 keys
        keys = list(set(keys))
        dfs = []
        prefix = ""
        for idx, key in enumerate(keys):
            resp = self.http.get(f"/experiment/{exp_id}/column/csv", params={"key": key})
            if resp.errmsg:
                continue

            url = resp.data.get("url", "")
            df = pd.read_csv(url, index_col=0)

            if idx == 0:
                # 从第一列名提取 prefix，例如 "t0707-02:17-loss_step" 中提取 "t0707-02:17-"
                first_col = df.columns[0]
                suffix = f"{key}_"
                if suffix in first_col:
                    prefix = first_col.split(suffix)[0]  # 结果为 "t0707-02:17-"
                else:
                    prefix = ""

            if prefix:
                df.columns = [
                    col[len(prefix):].removesuffix("_step") if col.startswith(prefix) else col.removesuffix("_step")
                    for col in df.columns
                ]
            else:
                df.columns = [col.removesuffix("_step") for col in df.columns]

            dfs.append(df)

        if not dfs:
            return ApiResponse[DataFrame](code=404, errmsg="No data found", data=pd.DataFrame())

        # 按列合并，使用 inner join 保证对齐 index
        result_df = pd.concat(dfs, axis=1, join="inner")

        return ApiResponse[DataFrame](code=200, errmsg="", data=result_df)
