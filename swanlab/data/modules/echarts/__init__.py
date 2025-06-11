"""
@author: ComPleHN
@file: __init__.py
@time: 2025/5/19 14:01
@desc: 集成 pyecharts
"""

import datetime

import pyecharts
import simplejson as json
from pyecharts.charts.base import Base
from pyecharts.commons import utils
from pyecharts.components import Table as T
from pyecharts.options import ComponentTitleOpts
from pyecharts.options.series_options import BasicOpts
from pyecharts.types import Sequence, Union, Optional
from swankit.core import MediaBuffer, DataSuite as D, MediaType

echarts = pyecharts.charts

PyEchartsBase = pyecharts.charts.base.Base
"""
pyecharts.charts.base.Base 的别名
"""

PyEchartsTable = T

__all__ = ["echarts", 'Echarts', 'PyEchartsBase', 'PyEchartsTable', "Table"]


class Table(T):
    def __init__(self):
        super().__init__()
        self.options: dict = {'_swanLab': "table", 'title_opts': None, 'headers': [], 'rows': []}

    def add(self, headers: Sequence, rows: Sequence, attributes: Optional[dict] = None):
        super().add(headers, rows, attributes)

        self.options.update({'headers': headers, 'rows': rows})
        return self

    def set_global_opts(self, title_opts: Union[ComponentTitleOpts, dict, None] = None):
        super().set_global_opts(title_opts)
        # Update title_opts in options
        if isinstance(title_opts, ComponentTitleOpts):
            self.options['title_opts'] = {
                'title': title_opts.title,
                'subtitle': title_opts.subtitle,
            }
        elif isinstance(title_opts, dict):
            self.options['title_opts'] = title_opts
        else:
            self.options['title_opts'] = None
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

    def dump_options(self) -> str:
        """序列化原始格式的options"""
        return utils.replace_placeholder(
            json.dumps(
                self.get_table_format(),
                default=self._default_parse,
                ignore_nan=True,
            )
        )


class Echarts(MediaType):
    def __init__(self, chart: Base):
        super().__init__()
        self._chart = chart
        self.buffer = MediaBuffer()

    # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        byte_string = self._chart.dump_options().encode('utf-8')
        hash_name = D.get_hash_by_bytes(byte_string)[:16]

        # 写入buffer
        self.buffer.write(byte_string)  # 写入二进制数据
        self.buffer.seek(0)  # 重置指针到开头，以便后续读取

        filename = f"echarts-step{self.step}-{hash_name}.json"
        return filename, self.buffer

    def get_chart(self):
        return self.Chart.ECHARTS

    def get_section(self):
        return "ECharts"
