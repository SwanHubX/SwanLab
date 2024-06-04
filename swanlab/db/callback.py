#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/12 15:22
@File: callback.py
@IDE: pycharm
@Description:
    注册回调函数
"""
from swanlab.data.run.callback import SwanLabRunCallback, ColumnInfo
from swanlab.log import swanlog
from .models import *
from .utils.chart import add_multi_chart
from .db_connect import connect
from .error import NotExistedError, ExistedError, ChartTypeError
from typing import Tuple, Callable, Optional
from swanlab.utils.file import check_exp_name_format, check_desc_format
from datetime import datetime
import time


class GlomCallback(SwanLabRunCallback):
    """
    GlomCallback类，swanlab本体与数据库的连接回调函数
    """

    def __init__(self):
        super(GlomCallback, self).__init__()
        self.exp: Optional[Experiment] = None

    @staticmethod
    def __get_exp_name(experiment_name: str = None, suffix: str = None) -> Tuple[str, str]:
        """
        预处理实验名称，如果实验名称过长，截断

        Parameters
        ----------
        experiment_name : str
            实验名称
        suffix : str
            实验名称后缀添加方式，可以为None、"default"或自由后缀，前者代表不添加后缀，后者代表添加时间戳后缀

        Returns
        ----------
        experiment_name : str
            校验后的实验名称
        exp_name : str
            最终的实验名称
        """
        # ---------------------------------- 校验实验名称 ----------------------------------
        experiment_name = "exp" if experiment_name is None else experiment_name
        # 校验实验名称
        experiment_name_checked = check_exp_name_format(experiment_name)
        # 如果实验名称太长的tip
        tip = "The experiment name you provided is too long, it has been truncated automatically."

        # 如果前后长度不一样，说明实验名称被截断了，提醒
        if len(experiment_name_checked) != len(experiment_name) and suffix is None:
            swanlog.warning(tip)

        # 如果suffix为None, 则不添加后缀，直接返回
        if suffix is None or suffix is False:
            return experiment_name_checked, experiment_name

        # 如果suffix为True, 则添加默认后缀
        if suffix is True:
            suffix = "default"

        # suffix必须是字符串
        if not isinstance(suffix, str):
            raise TypeError("The suffix must be a string, but got {}".format(type(suffix)))

        # 如果suffix_checked为default，则设置为默认后缀
        if suffix.lower().strip() == "default":
            # 添加默认后缀
            default_suffix = "{}".format(datetime.now().strftime("%b%d_%H-%M-%S"))
            exp_name = "{}_{}".format(experiment_name_checked, default_suffix)
        else:
            exp_name = "{}_{}".format(experiment_name_checked, suffix)

        # 校验实验名称，如果实验名称过长，截断
        experiment_name_checked = check_exp_name_format(exp_name)
        if len(experiment_name_checked) != len(exp_name):
            swanlog.warning(tip)

        return experiment_name_checked, exp_name

    def __str__(self):
        return "SwanLabGlomCallback"

    def on_init(self, proj_name: str, *args, **kwargs):
        # 连接本地数据库，要求路径必须存在，但是如果数据库文件不存在，会自动创建
        connect(autocreate=True)
        # 初始化项目数据库
        Project.init(proj_name)

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        suffix: str,
        setter: Callable[[str, str, str, str], None]
    ):
        # ---------------------------------- 实验描述校验 ----------------------------------
        description = description if description is not None else ""
        desc = check_desc_format(description)
        if desc != description:
            swanlog.warning("The description has been truncated automatically.")
        # ---------------------------------- 实验名称校验 ----------------------------------
        # 这个循环的目的是如果创建失败则等零点五秒重新生成后缀重新创建，直到创建成功
        # 但是由于需要考虑suffix为none不生成后缀的情况，所以需要在except中判断一下
        old = exp_name
        while True:
            experiment_name, exp_name = self.__get_exp_name(old, suffix)
            try:
                # 获得数据库实例
                self.exp = Experiment.create(name=exp_name, run_id=run_id, description=description, num=num)
                break
            except ExistedError:
                # 如果suffix名为default，说明是自动生成的后缀，需要重新生成后缀
                if isinstance(suffix, str) and suffix.lower().strip() == "default":
                    swanlog.debug(f"Experiment {exp_name} has existed, try another name...")
                    time.sleep(0.5)
                    continue
                # 其他情况下，说明是用户自定义的后缀，需要报错
                else:
                    Experiment.purely_delete(run_id=run_id)
                    raise ExistedError(f"Experiment {exp_name} has existed in local, please try another name.")
        # 执行相关设置
        setter(experiment_name, self.exp.light, self.exp.dark, self.exp.description)

    def on_log(self):
        # 每一次log的时候检查一下数据库中的实验状态
        # 如果实验状态不为0，说明实验已经结束，不允许再次调用log方法
        # 这意味着每次log都会进行查询，比较消耗性能，后续考虑采用多进程共享内存的方式进行优化
        swanlog.debug(f"Check experiment and state...")
        try:
            exp = Experiment.get(id=self.exp.id)
        except NotExistedError:
            raise KeyboardInterrupt("The experiment has been deleted by the user")
        # 此时self.__state == 0，说明是前端主动停止的
        if exp.status != 0:
            raise KeyboardInterrupt("The experiment has been stopped by the user")

    def on_column_create(self, column_info: ColumnInfo):
        # 创建Chart
        chart = Chart.create(
            column_info.key,
            experiment_id=self.exp,
            type=column_info.chart.value.chart_type,
            reference=column_info.reference,
            config=column_info.config,
        )
        # 创建命名空间，如果命名空间已经存在，会抛出ExistedError异常，捕获不处理即可
        # 需要指定sort，default命名空间的sort为0，其他命名空间的sort为None，表示默认添加到最后
        namespace = column_info.namespace
        try:
            n = Namespace.create(name=namespace, experiment_id=self.exp.id, sort=column_info.sort)
            swanlog.debug(f"Namespace {namespace} created, id: {n.id}")
        except ExistedError:
            n: Namespace = Namespace.get(name=namespace, experiment_id=self.exp.id)
            swanlog.debug(f"Namespace {namespace} exists, id: {n.id}")
        # 创建display，这个必然是成功的，因为display是唯一的，直接添加到最后一条即可
        Display.create(chart_id=chart.id, namespace_id=n.id)
        tag: Tag = Tag.create(
            experiment_id=self.exp.id,
            name=column_info.key,
            type=column_info.chart.value.chart_type,
        )
        # 添加一条source记录
        error = None
        if column_info.error is not None:
            error = {"data_class": column_info.error.got, "excepted": column_info.error.expected}
        Source.create(tag_id=tag.id, chart_id=chart.id, error=error)
        # 新建多实验对比图表数据
        try:
            add_multi_chart(tag_id=tag.id, chart_id=chart.id)
        except ChartTypeError:
            swanlog.warning("In the multi-experiment chart, the current type of tag is not as expected.")

    def on_stop(self, error: str = None):
        # 更新数据库中的实验状态
        self.exp.update_status(-1 if error is not None else 1)
