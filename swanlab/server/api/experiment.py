#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-03 00:39:18
@File: swanlab\server\api\experiment.py
@IDE: vscode
@Description:
    实验相关api，前缀：/experiment
"""
from fastapi import APIRouter
from ..utils import ResponseBody
from ...env import SWANLAB_LOGS_FOLDER
from ...database.project import ProjectTable
import os
import ujson

# from ...utils import create_time
from urllib.parse import unquote  # 转码路径参数
from typing import List, Dict
from ...utils import get_a_lock

router = APIRouter()

CONFIG_PATH = ProjectTable.path


# 获取当前实验信息
@router.get("/{experiment_id}")
async def get_experiment(experiment_id: int):
    """获取当前实验的信息

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    """
    # 读取 project.json 文件内容
    f = get_a_lock(CONFIG_PATH, "r")
    try:
        experiments: list = ujson.load(f)["experiments"]
        f.close()
        # 在experiments列表中查找对应实验的信息
        experiment = None
        for ex in experiments:
            if ex["experiment_id"] == experiment_id:
                experiment = ex
                break
        # 生成实验存储路径
        path = os.path.join(SWANLAB_LOGS_FOLDER, experiment["name"])
        experiment["tags"] = __list_subdirectories(path)
        return ResponseBody(0, data=experiment)
    except Exception as e:
        f.close()
        raise e


# # 修改实验的信息：名称/描述
# @router.patch("/{experiment_id}")
# async def patch_expriment(experiment_id: int, new_info: dict):
#     """修改当前实验的信息

#     parameter
#     ----------
#     experiment_id: int
#         实验唯一id，路径传参
#     new_info: dict
#         需要修改的内容
#         name: str
#         description: str
#     """
#     # 读取 project.json 文件内容
#     project: dict = ujson.load(open(CONFIG_PATH, "r", encoding="utf-8"))
#     experiments: list = project["experiments"]
#     # 在experiments列表中查找对应实验的信息
#     for index, experiment in enumerate(experiments):
#         if not experiment["experiment_id"] == experiment_id:
#             continue
#         if experiment["status"] == 0:
#             return ResponseBody(500, message="experiment not finished")
#         if new_info.get("name") is not None:
#             # 如果名字和之前一样
#             if new_info["name"] == experiment["name"]:
#                 return ResponseBody(500, message="experiment name is the same")
#             # 检查名字是否重复
#             if any(item["name"] == new_info["name"] for item in experiments):
#                 return ResponseBody(500, message="experiment name already exists")
#             # TODO 检查名字是否违规
#             # TODO 检查名字的步骤应该封装出来
#             # 修改实验目录的名字
#             os.rename(
#                 os.path.join(SWANLAB_LOGS_FOLDER, experiment["name"]),
#                 os.path.join(SWANLAB_LOGS_FOLDER, new_info["name"]),
#             )
#             # 在配置中修改
#             experiments[index]["name"] = new_info["name"]
#         if new_info.get("description") is not None:
#             experiment["description"] = new_info["description"]
#         project["experiments"] = experiments
#         project["update_time"] = create_time()
#         ujson.dump(project, open(CONFIG_PATH, "w", encoding="utf-8"), ensure_ascii=False, indent=4)
#         return ResponseBody(0)
#     return ResponseBody(500, message="experiment not found")


# 获取表单数据
@router.get("/{experiment_id}/{tag}")
async def get_tag_data(experiment_id: int, tag: str):
    """获取表单数据

    parameter
    ----------
    experiment_id: int
        实验唯一id，路径传参
    tag: str
        表单标签，路径传参，使用时需要 URIComponent 解码
    """
    tag = unquote(tag)
    # 读取实验信息内容
    experiments: list = ujson.load(open(CONFIG_PATH, "r", encoding="utf-8"))["experiments"]
    # 在experiments列表中查找对应实验的信息
    experiment_name = None
    for _, experiment in enumerate(experiments):
        if experiment["experiment_id"] == experiment_id:
            experiment_name = experiment["name"]
    # 找到最后一个还不存在，报错
    if experiment_name is None:
        # TODO 后续需要改成错误码
        raise KeyError(f'experiment id "{experiment_id}" not found')
    tag_data = __find_tag_data(experiment_name, tag)
    # 返回数据
    return ResponseBody(0, data={"sum": len(tag_data), "list": tag_data})


def __find_all_tag_data(base_path: str, paths: list) -> List[List[Dict]]:
    """读取path中所有的tag数据

    Parameters
    ----------
    base_path : str
        实验tag的存储路径
    paths : list
        tag的存储路径列表
    """
    tag_data_list: List[List[Dict]] = []
    for path in paths:
        # 读取tag数据，由于目前在设计上这些文件不会再被修改，所以不需要加锁
        with open(os.path.join(base_path, path), "r") as f:
            tag_data_list.append(ujson.load(f)["data"])
    return tag_data_list


def __find_tag_data(experiment_name: str, tag: str) -> (List[Dict], int):
    """找到对应实验的对应标签的数据

    Parameters
    ----------
    experiment_name : str
        实验名称
    tag : str
        tag的名称

    Returns
    -------
    List[Dict]
        tag的数据，内部是dict，包含每条tag的数据
        如果没有找到，返回空列表

    Raises
    ------
    KeyError
        如果tag不存在，抛出异常：'tag "{tag}" not found'
    """
    # 阈值，如果数据量大于阈值，只返回阈值条数据
    threshold = 5000
    # 获取tag对应的存储目录
    tag_path: str = os.path.join(SWANLAB_LOGS_FOLDER, experiment_name, tag)
    if not os.path.exists(tag_path):
        raise KeyError(f'tag "{tag}" not found')
    # 获取目录下存储的所有数据
    # 降序排列，最新的数据在最前面
    files: list = os.listdir(tag_path)
    if len(files) == 0:
        return []
    files.sort()
    tag_data: list = []
    # 最后一个文件代表当前数据量
    last_file = files[-1]
    tag_json = None
    # 锁住此文件，不再允许其他进程访问，换句话说，实验日志在log的时候被阻塞
    with get_a_lock(os.path.join(tag_path, last_file), mode="r") as f:
        # 读取数据
        tag_json = ujson.load(f)
        # 倒数第二个文件+当前文件的数据量等于总数据量
        # 倒数第二个文件可能不存在
        count = files[-2].split(".")[0] if len(files) > 1 else 0
        count = int(count) + len(tag_json["data"])
        # print(f"count={count}")
    # with生命周期结束，文件解锁，后续也不会再读取里面的数据了
    # 此时count代表总数据量，接下来按量倒叙读取数据
    if count <= threshold:
        # 读取所有数据
        # tag_json是最后一个文件的数据
        # 按顺序读取其他文件的数据
        tag_data_list = __find_all_tag_data(tag_path, files[:-1])
        # 将数据合并
        for data in tag_data_list:
            tag_data.extend(data)
        tag_data.extend(tag_json["data"])
        return tag_data
    else:
        # TODO 采样读取数据
        raise NotImplementedError("采样读取数据")


def __list_subdirectories(folder_path):
    # 使用 os.listdir 获取文件夹下所有项的列表
    items = os.listdir(folder_path)

    # 使用列表推导式筛选出所有的子文件夹
    subdirectories = [item for item in items if os.path.isdir(os.path.join(folder_path, item))]

    return subdirectories
