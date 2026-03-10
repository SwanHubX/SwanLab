"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API
"""


# def create_or_resume_experiment(
#     exp_name,
#     colors: Tuple[str, str],
#     description: Optional[str] = None,
#     job_type: Optional[str] = None,
#     group: Optional[str] = None,
#     tags: Optional[List[str]] = None,
#     created_at: Optional[str] = None,
#     cuid: Optional[str] = None,
#     must_exist: bool = False,
# ) -> bool:
#     """
#     初始化实验，获取存储信息
#     :param exp_name: 所属实验名称
#     :param colors: 实验颜色，有两个颜色
#     :param description: 实验描述
#     :param job_type: 任务类型
#     :param group: 实验组
#     :param tags: 实验标签
#     :param created_at: 实验创建时间，格式为 ISO 8601
#     :param cuid: 实验的唯一标识符，如果不提供则由后端生成
#     :param must_exist: 如果 cuid 被传递，是否限制实验必须存在
#
#     :raises RuntimeError: 如果实验不存在且must_exist为True
#     :raises NotImplementedError: 如果项目未挂载
#
#     :return: 返回实验为新建的还是更新的，为 True 时为新建实验
#     """
#     if must_exist:
#         assert cuid is not None, "cuid must be provided when must_exist is True"
#         try:
#             self.get(f"/project/{self.groupname}/{self.__proj.name}/runs/{cuid}")
#         except ApiError as e:
#             if e.resp.status_code == 404 and e.resp.reason == "Not Found":
#                 raise RuntimeError(f"Experiment {cuid} does not exist in project {self.projname}")
#
#     labels = [{"name": tag} for tag in tags] if tags else []
#     post_data = {
#         "name": exp_name,
#         "description": description,
#         "createdAt": created_at,
#         "colors": list(colors),
#         "labels": labels if len(labels) else None,
#         "job": job_type,
#         "cluster": group,
#         "cuid": cuid,
#     }
#     post_data = {k: v for k, v in post_data.items() if v is not None}  # 移除值为None的键
#
#     # 这部分错误将不会被上层捕获，直接抛出异常
#     try:
#         data, resp = self.post(f"/project/{self.groupname}/{self.__proj.name}/experiment", post_data)
#     except ApiError as e:
#         if e.resp.status_code == 400 and e.resp.reason == "Bad Request":
#             # 指定的 cuid 对应的实验是克隆实验
#             raise ValueError(
#                 f"Experiment with CUID {cuid} is a cloned experiment (cloned experiments cannot be resumed).",
#             )
#         elif e.resp.status_code == 403 and e.resp.reason == "Forbidden":
#             # 权限不足
#             raise ValueError(f"Project permission denied: {self.projname}")
#         elif e.resp.status_code == 404 and e.resp.reason == "Not Found":
#             # 传入的项目不存在
#             raise ValueError(f"Project {self.projname} not found")
#         elif e.resp.status_code == 404 and e.resp.reason == "Disabled Resource":
#             # 传入的实验被删除
#             raise ValueError(f"Experiment {cuid} has been deleted")
#         elif e.resp.status_code == 409 and e.resp.reason == "Conflict":
#             # 传入 cuid 但是实验不属于当前项目
#             raise ValueError(f"Experiment with CUID {cuid} does not belong to project {self.projname}")
#         raise e
#     return new
