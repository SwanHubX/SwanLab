#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-01-01 15:58:58
@File: swanlab/data/run/main.py
@IDE: vscode
@Description:
    在此处定义SwanLabRun类并导出
"""
from ..settings import SwanDataSettings, get_runtime_project
from ...log import register, swanlog
from ..system import get_system_info
from .utils import (
    get_a_lock,
    check_exp_name_format,
    check_desc_format,
    get_package_version,
    create_time,
    generate_color,
)
from datetime import datetime
import sys, os
import random
import ujson
from .exp import SwanLabExp


class SwanLabRun:
    """
    The SwanLabRun class is used for logging during a single experiment.
    There should be only one instance of the SwanLabRun class for each experiment.
    """

    def __init__(
        self,
        experiment_name: str = None,
        description: str = None,
        config: dict = None,
        log_level: str = None,
        suffix: str = None,
    ):
        """
        Initializing the SwanLabRun class involves configuring the settings and initiating other logging processes.

        Parameters
        ----------
        experiment_name : str, optional
            实验名称，实验名称应该唯一，由0-9，a-z，A-Z，" ","_","-","/"组成
            如果不提供此参数(为None)，SwanLab将自动生成一个实验名称
        description : str, optional
            实验描述，用于对当前实验进行更详细的介绍或标注
            如果不提供此参数(为None)，可以在web界面中进行修改,这意味着必须在此改为空字符串""
        config : dict, optional
            实验参数配置，可以在web界面中显示，如学习率、batch size等
            不需要做任何限制，但必须是字典类型，可被json序列化，否则会报错
        log_level : str, optional
            当前实验的日志等级，默认为'info'，可以从'debug'、'info'、'warning'、'error'、'critical'中选择
            不区分大小写，如果不提供此参数(为None)，则默认为'info'
            如果提供的日志等级不在上述范围内，默认改为info
        suffix : str, optional
            实验名称后缀，用于区分同名实验，格式为yyyy-mm-dd_HH-MM-SS
            如果不提供此参数(为None)，不会添加后缀
        """
        # ---------------------------------- 初始化类内参数 ----------------------------------
        # 生成一个唯一的id，随机生成一个8位的16进制字符串，小写
        self.__id = hex(random.randint(0, 2**32 - 1))[2:].zfill(8)
        self.__timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        # ---------------------------------- 初始化实验配置 ----------------------------------
        # FIXME 初始化实验配置这里有点问题，写入的时候分了两个步骤
        # 等后续使用数据库后，SwanDataSettings传入的应该是当前run对象的str表示，而不是实验名称
        # 确保project.json文件存在
        project_path = get_runtime_project()
        open(project_path, "a", encoding="utf-8").close()
        exp_name, cut = self.__get_exp_name(experiment_name, project_path, suffix)
        # ---------------------------------- 初始化日志等级 ----------------------------------
        level = self.__check_log_level(log_level)
        # ---------------------------------- 初始化其他对象 ----------------------------------
        # 初始化配置
        self.settings = SwanDataSettings(exp_name)
        # 初始化日志记录器
        # FIXME 日志记录器的初始化强依赖于settings，这不行，等后续使用数据库后，在data部分的顶层完成日志记录器的初始化
        register(self.settings.output_path, self.settings.console_dir, log_level=level)
        # 在此处判断是否被截断了，warning一下
        if cut:
            swanlog.warning(f"The experiment name you provided is too long, it has been truncated to {exp_name}.")
        # 并判断一下log_level是否合法，如果不合法，warning一下
        if log_level != level and log_level is not None:
            swanlog.warning(f"The log level you provided is not valid, it has been set to {level}.")
        # ---------------------------------- 其余配置 ----------------------------------
        # 注册实验
        self.exp = self.__register_exp(
            exp_name,
            self.__check_description(description),
            self.__check_config(config),
        )

    def __str__(self) -> str:
        """此类的字符串表示"""
        return "run-{}-{}".format(self.__timestamp, self.__id)

    def log(self, data: dict, step: int = None):
        if not isinstance(data, dict):
            return swanlog.error(
                "log data must be a dict, but got {}, SwanLab will ignore records it.".format(type(data))
            )
        # 检查step的类型
        if step is not None and (not isinstance(step, int) or step < 0):
            swanlog.error(
                "'step' must be an integer not less than zero, but got {}, SwanLab will automatically set step".format(
                    step
                )
            )
            step = None
        # 遍历data，记录data
        for key in data:
            # 遍历字典的key，记录到本地文件中
            d = data[key]
            # 数据类型的检查将在创建chart配置的时候完成，因为数据类型错误并不会影响实验进行
            self.exp.add(key, d, step=step)

    def success(self):
        """标记实验成功"""
        self.__set_exp_status(1)

    def fail(self):
        """标记实验失败"""
        self.__set_exp_status(-1)

    def __set_exp_status(self, status: int):
        if status not in [1, -1, 0]:
            swanlog.warning("Invalid status when set, status must be 1, -1 or 0, but got {}".format(status))
            swanlog.warning("SwanLab will set status to -1")
            status = -1
        # 锁上文件，更新实验状态
        with get_a_lock(self.settings.project_path, "r+") as file:
            project = ujson.load(file)
            for index, experiment in enumerate(project["experiments"]):
                if experiment["experiment_id"] == self.exp.id:
                    project["experiments"][index]["status"] = status
                    project["update_time"] = create_time()
                    break
            file.truncate(0)
            file.seek(0)
            ujson.dump(project, file)

    def __get_exp_name(self, experiment_name: str = None, project_path: str = None, suffix: str = None) -> tuple:
        """拿到实验名称

        Parameters
        ----------
        experiment_name : str
            实验名称
        cut : bool
            是否被截断
        """
        max_len = 20
        cut = experiment_name is not None and len(experiment_name) > max_len
        experiment_name = "exp" if experiment_name is None else check_exp_name_format(experiment_name)
        # 为实验名称添加后缀，格式为yyyy-mm-dd_HH-MM-SS
        if suffix is not None and suffix.lower() != "timestamp":
            suffix = "timestamp"
            swanlog.warning(f"The suffix you provided is not valid, it has been set to {suffix}.")
        # ---------------------------------- 自动添加后缀 ----------------------------------
        if suffix is not None:
            suffix = self.__timestamp
            experiment_name = f"{experiment_name}_{suffix}"
        # 拿取project.json文件中的实验配置
        with get_a_lock(project_path, "r+") as f:
            # 如果project.json文件为空，创建一个新的project.json文件
            project_exist = os.path.exists(project_path) and os.path.getsize(project_path) != 0
            if not project_exist:
                return experiment_name, cut
            project = ujson.load(f)
            # 检查实验名称是否存在
            while True:
                unique = True
                # 遍历所有实验，检查实验名称是否存在
                for exp in project["experiments"]:
                    if exp["name"] == experiment_name:
                        # 第一次遇到重名，添加后缀-1
                        if unique:
                            experiment_name += "-1"
                            unique = False
                        # 第二次遇到重名，后缀+1
                        else:
                            num = int(experiment_name.split("-")[-1]) + 1
                            experiment_name = "".join(experiment_name.split("-")[:-1]) + f"-{num}"
                            unique = False
                        break
                # 如果实验名称不存在，跳出循环
                if unique:
                    break
                swanlog.warning(
                    f"The experiment name you provided is not unique, it has been set to {experiment_name}, check again."
                )
        return experiment_name, cut

    def __register_exp(
        self,
        experiment_name: str,
        description: str = None,
        config: dict = None,
    ) -> SwanLabExp:
        """将此实验配置写入project.json文件中
        如果project.json文件不存在，则创建一个
        最终返回实验名称
        """
        project_path = self.settings.project_path
        with get_a_lock(project_path, "r+") as f:
            # 如果project.json文件为空，创建一个新的project.json文件
            project_exist = os.path.exists(project_path) and os.path.getsize(project_path) != 0
            project = ujson.load(f) if project_exist else self.__new_project()
            # 创建新的实验配置
            project["_sum"] = project["_sum"] + 1
            experiment = self.__new_experiment(project["_sum"], experiment_name, description, config)
            # 添加实验配置到project.json文件中
            project["experiments"].append(experiment)
            f.truncate(0)
            f.seek(0)
            ujson.dump(project, f)
        swanlog.info(f"Experiment {experiment_name} has been registered.")
        return SwanLabExp(self.settings, experiment["experiment_id"])

    def __check_log_level(self, log_level: str) -> str:
        """检查日志等级是否合法"""
        valids = ["debug", "info", "warning", "error", "critical"]
        if log_level is None:
            return "info"
        elif log_level.lower() in valids:
            return log_level.lower()
        else:
            return "info"

    def __check_description(self, description: str) -> str:
        """检查实验描述是否合法"""
        if description is None:
            return ""
        desc = check_desc_format(description)
        if desc != description:
            swanlog.warning("The description has been truncated automatically.")
        return desc

    def __check_config(self, config: dict) -> dict:
        """检查实验配置是否合法"""
        if config is None:
            return {}
        # config必须可以被json序列化
        try:
            ujson.dumps(config)
        except:
            raise TypeError(f"config: {config} is not a valid dict, which can be json serialized")
        return config

    def __new_project(self):
        """创建一个新的project.json文件"""
        time = create_time()
        return {
            "name": os.path.basename(self.settings.root_dir),
            "_sum": 0,
            "experiments": [],
            "version": get_package_version(),
            "create_time": time,
            "update_time": time,
        }

    @staticmethod
    def __new_experiment(sum: int, name: str, description: str, config: dict):
        """创建一个新的实验配置"""
        return {
            "experiment_id": sum,
            "name": name,
            "version": get_package_version(),
            "status": 0,
            "description": description,
            "config": config,
            "color": generate_color(sum),
            "system": get_system_info(),
            "argv": sys.argv,
            "index": sum,
            "create_time": create_time(),
            "update_time": create_time(),
        }
