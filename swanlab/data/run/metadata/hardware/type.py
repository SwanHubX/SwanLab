"""
@author: cunyue
@file: type.py
@time: 2024/12/3 20:27
@description: 硬件信息采集类型定义
"""

from abc import ABC, abstractmethod
from typing import TypedDict, Tuple, Optional, Any, List, Union

from swankit.callback.models import ColumnConfig, YRange

from swanlab.data.run.namer import generate_colors
from swanlab.log import swanlog


# 定义硬件信息类型
class HardwareInfo(TypedDict):
    # 硬件信息键
    key: str
    # 硬件信息值
    value: Union[str, int, float]
    # 硬件信息名称
    name: str
    # 相关配置
    config: Optional[ColumnConfig]


HardwareInfoList = List[HardwareInfo]


class HardwareConfig(ColumnConfig):
    """
    继承自ColumnConfig的硬件配置类
    每clone一次，metric_color都会自动生成一个
    如果不指定，初次生成的metric_color为None
    """

    def __init__(self, y_range=None, chart_name=None, chart_index=None, metric_name=None, metric_color=None):
        self.__cloned = 0
        super().__init__(
            y_range=y_range,
            chart_name=chart_name,
            chart_index=chart_index,
            metric_name=metric_name,
            metric_color=metric_color,
        )

    def clone(
        self,
        y_range: YRange = None,
        chart_name: Optional[str] = None,
        chart_index: Optional[str] = None,
        metric_name: Optional[str] = None,
        metric_color: Optional[Tuple[str, str]] = None,
    ):
        """
        重写clone方法，每次clone都会生成一个新的metric_color
        :param y_range: y轴范围
        :param chart_name: 图表名称
        :param chart_index: 图表索引
        :param metric_name: 指标名称
        :param metric_color: 指标颜色
        :return: 新的HardwareConfig对象
        """
        try:
            if metric_color is None:
                metric_color = generate_colors(self.__cloned)
            return HardwareConfig(
                y_range=y_range if y_range is not None else self.y_range,
                chart_name=chart_name if chart_name is not None else self.chart_name,
                metric_name=metric_name if metric_name is not None else self.metric_name,
                chart_index=chart_index if chart_index is not None else self.chart_index,
                metric_color=metric_color,
            )
        finally:
            self.__cloned += 1


class CollectGuard:
    """
    采集守卫，在采集任务执行前后可以选择执行一些操作
    这些操作可能比较耗费资源，因此需要选择性执行，而选择的依据是是执行的次数
    这与任务的采集间隔有关，即此类适配了 MonitorCron 的采集间隔
    """

    def __init__(self):
        self.collect_num = 0

    def before_collect(self):
        try:
            # count=0，执行一次
            if self.collect_num == 0:
                return self.before_collect_impl()
            # 60>count>0, 不执行
            elif self.collect_num < 60:
                return
            else:
                return self.before_collect_impl()
        finally:
            self.collect_num += 1

    def after_collect(self):
        # 60>count>0, 不执行
        if self.collect_num < 60:
            return
        # count=60，执行一次
        elif self.collect_num == 60:
            return self.after_collect_impl()
        else:
            return self.after_collect_impl()

    def before_collect_impl(self):
        pass

    def after_collect_impl(self):
        pass


class HardwareCollector(CollectGuard, ABC):
    @abstractmethod
    def collect(self) -> HardwareInfoList:
        """
        采集硬件信息的主函数，子类覆写此方法实现具体的硬件信息采集业务逻辑
        """
        pass

    @staticmethod
    def try_run():
        """
        一个装饰器
        可以选择对某个采集任务的单个函数进行异常捕获，避免因为某个采集任务失败导致整个采集任务失败
        """

        def decorator(func):
            def wrapper(*args, **kwargs):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    swanlog.debug(f"Atoms collection failed: {func.__name__}, {str(e)}")
                    return None

            return wrapper

        return decorator

    def __call__(self) -> Optional[HardwareInfoList]:
        """
        聚合执行采集任务，包括采集任务的前置和后置操作
        """
        try:
            self.before_collect()
            result = self.collect()
            if result is None:
                return None
            # 过滤掉采集结果中的None
            return [r for r in result if r is not None]
        except NotImplementedError as n:
            raise n
        except Exception as e:
            swanlog.debug(f"Collection failed: {self.__class__.__name__}, {str(e)}")
            return None
        finally:
            self.after_collect()


# 定义硬件信息执行函数的返回结果
HardwareFuncResult = Tuple[Optional[Any], Optional[HardwareCollector]]
