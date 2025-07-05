from typing import Dict, Optional

from swanlab.data.modules import DataWrapper
from swanlab.log import swanlog
from swanlab.toolkit import (
    MetricInfo,
    MetricErrorInfo,
    ColumnClass,
    SectionType,
    ColumnConfig,
)
from .helper import SwanLabRunOperator
from .key import SwanLabKey
from ..store import get_run_store


class SwanLabExp:
    """
    Class for running experiments
    save keys when running experiments
    """

    def __init__(self, operator: SwanLabRunOperator) -> None:
        """初始化实验

        Parameters
        ----------
        operator : SwanLabRunOperator
            操作员
        """
        self._run_store = get_run_store()
        self._operator = operator
        # 当前实验的所有tag数据字段
        self._keys: Dict[str, SwanLabKey] = {}
        # 恢复实验时，同步云端实验的指标数据
        if self._run_store.metrics is not None:
            media_dir, log_dir = self._run_store.media_dir, self._run_store.log_dir
            for kid, (key, (column_type, column_class, error, step)) in enumerate(self._run_store.metrics.items()):
                key_index = self._generate_key_index(key, column_class)
                self._keys[key_index], column_info = SwanLabKey.mock_from_remote(
                    key, column_type, column_class, error, media_dir, log_dir, kid, step
                )
                self._operator.on_column_create(column_info)

    @staticmethod
    def _generate_key_index(key: str, column_class: ColumnClass) -> str:
        """
        生成key的索引
        :param key: key的云端标识
        :param column_class: 列类型，CUSTOM为自定义key，SYSTEM为系统key
        """
        if column_class not in ('CUSTOM', 'SYSTEM'):
            raise RuntimeError(
                f"Column class must be 'CUSTOM' or 'SYSTEM', got {column_class}. "
                f"Maybe you need to update swanlab: pip install -U swanlab"
            )
        return f"{column_class}-{key}"

    def _warn_type_error(self, key_index: str, key: str):
        """警告类型错误
        执行此方法时需保证key已经存在
        """
        key_obj = self._keys[key_index]
        class_name = key_obj.column_info.got
        expected = key_obj.column_info.expected
        if key_obj.is_chart_valid:
            return
        if class_name == "list":
            swanlog.error(f"Data type error, key: {key}, there is element of invalid data type in the list.")
        else:
            swanlog.error(f"Data type error, key: {key}, data type: {class_name}, expected: {expected}")

    def _warn_chart_error(self, key_index: str, key: str):
        """
        警告图表创建错误
        执行此方法时需保证key已经存在
        """
        key_obj = self._keys[key_index]
        if key_obj.is_chart_valid:
            return
        class_name = key_obj.column_info.got
        if class_name == "list":
            swanlog.warning(
                f"Chart '{key}' creation failed. "
                f"Reason: The data type in list of the key '{key}' is not as expected, please check the data type."
            )
        else:
            swanlog.warning(
                f"Chart '{key}' creation failed. "
                f"Reason: The expected value type for the chart '{key}' is one of int,"
                f"float or BaseType, but the input type is {class_name}."
            )

    def _add(
        self,
        key: str,
        name: Optional[str],
        column_class: ColumnClass,
        column_config: Optional[ColumnConfig],
        section_type: SectionType,
        data: DataWrapper,
        step: int = None,
    ) -> MetricInfo:
        """记录一条新的key数据

        Parameters
        ----------
        key : str
            key的云端唯一标识
        name : str
            key的实际名称
        column_class : str
            列类型，CUSTOM为自定义key，SYSTEM为系统key
        section_type : str
            key的组类型
        data : DataWrapper
            包装后的数据，用于数据解析
        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        key_index = self._generate_key_index(key, column_class)
        # 判断tag是否存在，如果不存在则创建tag
        key_obj: SwanLabKey = self._keys.get(key_index, None)

        # ---------------------------------- 包装器解析 ----------------------------------

        if step is not None and not isinstance(step, int):
            swanlog.warning(f"Step {step} is not int, SwanLab will set it automatically.")
            step = None
        if key_obj is None:
            step = 0 if step is None or not isinstance(step, int) else step
        else:
            step = len(key_obj.steps) if step is None else step
            if step in key_obj.steps:
                swanlog.debug(f"Step {step} on key {key} already exists, ignored.")
                return MetricErrorInfo(column_info=key_obj.column_info, error=DataWrapper.create_duplicate_error())
        data.parse(step=step, key=key)

        # ---------------------------------- 图表创建 ----------------------------------

        if key_obj is None:
            num = len(self._keys)
            # 将此tag对象添加到实验列表中
            key_obj = SwanLabKey(key, self._run_store.media_dir, self._run_store.log_dir)
            self._keys[key_index] = key_obj
            # 新建图表，完成数据格式校验
            column_info = key_obj.create_column(
                key,
                name,
                column_class,
                column_config,
                section_type,
                data,
                num,
            )
            self._warn_type_error(key_index, key)
            # 创建新列，生成回调
            self._operator.on_column_create(column_info)

        # 检查tag创建时图表是否创建成功，如果失败则也没有写入数据的必要了，直接退出
        if not key_obj.is_chart_valid:
            self._warn_chart_error(key_index, key)
            return MetricErrorInfo(key_obj.column_info, error=key_obj.column_info.error)
        key_info = key_obj.add(data)
        key_info.buffers = data.parse().buffers
        key_info.media_dir = self._run_store.media_dir
        return key_info

    def add(
        self,
        data: DataWrapper,
        key: str,
        name: str = None,
        column_class: ColumnClass = 'CUSTOM',
        column_config: Optional[ColumnConfig] = None,
        section_type: SectionType = "PUBLIC",
        step: int = None,
    ) -> MetricInfo:
        """记录一条新的key数据
        Parameters
        ----------
        data : DataWrapper
            包装后的数据，用于数据解析
        key : str
            列的云端唯一标识
        name : str
            列的实际名称, 默认与key相同
        column_class : str, optional
            列的类型
        column_config : Optional[ColumnConfig], optional
            列的额外配置信息
        section_type : str, optional
            key的组类型
        step : int, optional
            步数，如果不传则默认当前步数为'已添加数据数量+1'
            在log函数中已经做了处理，此处不需要考虑数值类型等情况
        """
        section_type: SectionType = 'CUSTOM' if data.is_custom else section_type
        m = self._add(key, name, column_class, column_config, section_type, data, step)
        self._operator.on_metric_create(m)
        return m
