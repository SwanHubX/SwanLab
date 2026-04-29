"""
@author: cunyue
@file: __init__.py.py
@time: 2026/4/29 18:46
@description: swanlab.echarts 子模块，导出 pyecharts 图表类型

用户可通过 swanlab.echarts.Bar()、swanlab.echarts.Line() 等创建图表，
然后用 swanlab.log({"key": chart}) 记录。
"""

from pyecharts import options
from pyecharts.charts import (
    Bar,
    Bar3D,
    Boxplot,
    Calendar,
    EffectScatter,
    Funnel,
    Gauge,
    Graph,
    Grid,
    HeatMap,
    Kline,
    Line,
    Line3D,
    Liquid,
    Parallel,
    PictorialBar,
    Pie,
    Polar,
    Radar,
    Sankey,
    Scatter,
    Scatter3D,
    Sunburst,
    Surface3D,
    ThemeRiver,
    Tree,
    TreeMap,
)
from pyecharts.charts.base import Base

from .metrics import confusion_matrix, pr_curve, roc_curve
from .table import Table

__all__ = [
    "options",
    "Base",
    "Bar",
    "Bar3D",
    "Boxplot",
    "Calendar",
    "EffectScatter",
    "Funnel",
    "Gauge",
    "Graph",
    "Grid",
    "HeatMap",
    "Kline",
    "Line",
    "Line3D",
    "Liquid",
    "Parallel",
    "PictorialBar",
    "Pie",
    "Polar",
    "Radar",
    "Sankey",
    "Scatter",
    "Scatter3D",
    "Sunburst",
    "Surface3D",
    "ThemeRiver",
    "Tree",
    "TreeMap",
    "Table",
    "roc_curve",
    "pr_curve",
    "confusion_matrix",
]
