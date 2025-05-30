import json

import swanlab
from swanlab.data.modules import Echarts


def test_echarts():
    """测试图表功能"""
    chart = swanlab.echarts.Bar().add_xaxis(["X", "Y"]).add_yaxis("数值", [10, 20])
    assert isinstance(chart, swanlab.echarts.Bar)
    c = Echarts(chart)
    assert isinstance(c, Echarts)
    filename, buffer = c.parse()
    # 返回文件名称
    assert isinstance(filename, str)
    assert filename.endswith(".json")
    assert buffer is not None
    # 可解析为字典
    chart_data = json.loads(buffer.getvalue().decode("utf-8"))
    assert isinstance(chart_data, dict)
