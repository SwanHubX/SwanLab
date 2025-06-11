"""
@author: cunyue
@file: table.py
@time: 2025/6/11 11:43
@description: 额外适配 pyecharts 的 table 类型，前端重写，自定义数据结构
"""

import datetime

import simplejson as json
from pyecharts.commons import utils
from pyecharts.components import Table as T
from pyecharts.options import ComponentTitleOpts
from pyecharts.options.series_options import BasicOpts
from pyecharts.types import Sequence, Union, Optional


class Table(T):
    def __init__(self):
        super().__init__()
        self.options: dict = {'_swanLab': "table", 'headers': [], 'rows': []}

    def add(self, headers: Sequence, rows: Sequence, attributes: Optional[dict] = None):
        super().add(headers, rows, attributes)
        self.options.update({'headers': headers, 'rows': rows})
        return self

    def get_table_format(self) -> dict:
        """获取转换后的表格格式(包含rowData和colDefs)"""
        original_data = utils.remove_key_with_none_value(self.options)

        # 创建新字典，保留所有原始字段
        converted_data = original_data.copy()

        # 转换表格数据部分
        if "headers" in original_data and "rows" in original_data:
            # 转换列定义
            converted_data["colDefs"] = [{"field": header} for header in original_data["headers"]]

            # 转换行数据
            converted_data["rowData"] = [dict(zip(original_data["headers"], row)) for row in original_data["rows"]]

            # 移除原始的headers和rows字段
            converted_data.pop("headers", None)
            converted_data.pop("rows", None)

        return converted_data

    def set_global_opts(self, title_opts: Union[ComponentTitleOpts, dict, None] = None):
        raise NotImplementedError("set_global_opts is not supported in swanlab.echarts.Table")

    @staticmethod
    def _default_parse(o):
        """
        默认的序列化方法，处理日期、JsCode和BasicOpts等特殊类型
        """
        if isinstance(o, (datetime.date, datetime.datetime)):
            return o.isoformat()
        if isinstance(o, utils.JsCode):
            return o.replace("\\n|\\t", "").replace(r"\\n", "\n").replace(r"\\t", "\t").js_code
        if isinstance(o, BasicOpts):
            if isinstance(o.opts, Sequence):
                return [utils.remove_key_with_none_value(item) for item in o.opts]
            else:
                return utils.remove_key_with_none_value(o.opts)
        return str(o)  # 对于其他类型，转换为字符串

    def dump_options(self) -> str:
        """序列化原始格式的options"""
        return utils.replace_placeholder(
            json.dumps(self.get_table_format(), default=self._default_parse, ignore_nan=True)
        )
