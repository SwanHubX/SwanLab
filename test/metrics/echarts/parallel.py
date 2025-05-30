"""
// @author: ComPleHN
// @file: parallel.vue
// @time: 2025/5/28 12:54
// @description: 本文件是对于echarts的 平行坐标系 图表测试
"""
# ---------------------------------------------- Parallel - Basic_parallel ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Parallel

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=parallel-simple

目前无法实现的功能:

1、
"""

parallel_axis = [
    {"dim": 0, "name": "Price"},
    {"dim": 1, "name": "Net Weight"},
    {"dim": 2, "name": "Amount"},
    {
        "dim": 3,
        "name": "Score",
        "type": "category",
        "data": ["Excellent", "Good", "OK", "Bad"],
    },
]

data = [[12.99, 100, 82, "Good"], [9.99, 80, 77, "OK"], [20, 120, 60, "Excellent"]]


c1 = (
    Parallel()
    .add_schema(schema=parallel_axis)
    .add(
        series_name="",
        data=data,
        linestyle_opts=opts.LineStyleOpts(width=4, opacity=0.5),
    )
)

# ---------------------------------------------- Parallel - Parallel_category ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Parallel

data = [
    [1, 91, 45, 125, 0.82, 34, 23, "良"],
    [2, 65, 27, 78, 0.86, 45, 29, "良"],
    [3, 83, 60, 84, 1.09, 73, 27, "良"],
    [4, 109, 81, 121, 1.28, 68, 51, "轻度污染"],
    [5, 106, 77, 114, 1.07, 55, 51, "轻度污染"],
    [6, 109, 81, 121, 1.28, 68, 51, "轻度污染"],
    [7, 106, 77, 114, 1.07, 55, 51, "轻度污染"],
    [8, 89, 65, 78, 0.86, 51, 26, "良"],
    [9, 53, 33, 47, 0.64, 50, 17, "良"],
    [10, 80, 55, 80, 1.01, 75, 24, "良"],
    [11, 117, 81, 124, 1.03, 45, 24, "轻度污染"],
    [12, 99, 71, 142, 1.1, 62, 42, "良"],
    [13, 95, 69, 130, 1.28, 74, 50, "良"],
    [14, 116, 87, 131, 1.47, 84, 40, "轻度污染"],
]
c2 = (
    Parallel()
    .add_schema(
        [
            opts.ParallelAxisOpts(dim=0, name="data"),
            opts.ParallelAxisOpts(dim=1, name="AQI"),
            opts.ParallelAxisOpts(dim=2, name="PM2.5"),
            opts.ParallelAxisOpts(dim=3, name="PM10"),
            opts.ParallelAxisOpts(dim=4, name="CO"),
            opts.ParallelAxisOpts(dim=5, name="NO2"),
            opts.ParallelAxisOpts(dim=6, name="CO2"),
            opts.ParallelAxisOpts(
                dim=7,
                name="等级",
                type_="category",
                data=["优", "良", "轻度污染", "中度污染", "重度污染", "严重污染"],
            ),
        ]
    )
    .add("parallel", data)
    .set_global_opts(title_opts=opts.TitleOpts(title="Parallel-Category"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="parallel",
    public=True,
)

swanlab.log(
    {
        "Parallel - Basic_parallel": c1,
        "Parallel - Parallel_category": c2
    }
)
