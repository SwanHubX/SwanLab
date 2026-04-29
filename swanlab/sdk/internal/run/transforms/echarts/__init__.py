"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 19:17
@description: ECharts TransformMedia 子类，将 pyecharts 图表转换为 JSON 文件并存储
"""

import hashlib
from pathlib import Path

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.echarts import EChartsDataType


class ECharts(TransformMedia):
    def __init__(self, chart: EChartsDataType, caption: CaptionType = None):
        """
        Parameters
        ----------
        chart : pyecharts 图表对象
            支持 pyecharts.charts.* 的所有图表类型以及 pyecharts.components.Table。
            对象必须实现 dump_options() 方法。
        caption : str, optional
            图表说明
        """
        super().__init__()
        attrs = self._unwrap(chart)
        if attrs:
            self.json_content: str = attrs["json_content"]
            self.caption = caption if caption is not None else attrs.get("caption")
            return

        if not callable(getattr(chart, "dump_options", None)):
            raise TypeError(
                f"Unsupported chart type: {type(chart).__name__}. "
                "Expected a pyecharts chart object with dump_options() method."
            )
        self.json_content = chart.dump_options()  # type: ignore
        self.caption = caption

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_ECHARTS

    def transform(self, *, step: int, path: Path) -> MediaItem:
        content_encode = self.json_content.encode("utf-8")
        sha256 = hashlib.sha256(content_encode).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.json"
        # safe_write 默认 mode="w"（文本模式），Windows 上会将 \n 转换为 \r\n。而 ECharts.transform 中
        # json_content.encode("utf-8") 计算的 size 和 sha256 是基于 \n 的。
        fs.safe_write(path / filename, content_encode, mode="wb")
        return MediaItem(filename=filename, sha256=sha256, size=len(content_encode), caption=self.caption or "")
