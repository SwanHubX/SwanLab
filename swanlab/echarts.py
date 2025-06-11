"""
@author: cunyue
@file: swanlab_echarts.py
@time: 2025/6/11 17:27
@description: echarts 类型是 swanlab 数据类的一部分，为了方便用户使用 echarts 图表，在 swanlab 包顶层新建 echarts 模块
使 echarts 支持如下导出模式：
```python
from swanlab.echarts import Table
from swanlab.echarts import option as eopt
```

反之，如果不建此文件，那么就只能：
```python
from swanlab import echarts

Table = echarts.Table
eopt = echarts.option
```

此外，由于前端支持有限，这里仅导出部分 echarts 类型
"""

from pyecharts import options
from pyecharts.charts import *

from .data.modules.custom_charts.table import Table

__all__ = [
    "options",
    "Bar3D",
    "Bar",
    "Boxplot",
    "Calendar",
    "Kline",
    "Grid",
    "Scatter",
    "EffectScatter",
    "Funnel",
    "Gauge",
    "Graph",
    "HeatMap",
    "Line",
    "Line3D",
    "Liquid",
    "Parallel",
    "PictorialBar",
    "Pie",
    "Polar",
    "Radar",
    "Sankey",
    "Scatter3D",
    "Sunburst",
    "Surface3D",
    "ThemeRiver",
    "Tree",
    "TreeMap",
    "Table",
]
