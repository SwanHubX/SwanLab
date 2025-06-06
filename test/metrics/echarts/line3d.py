"""
// @author: ComPleHN
// @file: line3d.vue
// @time: 2025/5/27 18:39
// @description: 本文件是对于echarts的 3D折线图 图表测试 
"""
# ---------------------------------------------- Line3d - Line3d_autorotate ----------------------------------------------
import math

from pyecharts import options as opts
from pyecharts.charts import Line3D
from pyecharts.faker import Faker

data = []
for t in range(0, 25000):
    _t = t / 1000
    x = (1 + 0.25 * math.cos(75 * _t)) * math.cos(_t)
    y = (1 + 0.25 * math.cos(75 * _t)) * math.sin(_t)
    z = _t + 2.0 * math.sin(75 * _t)
    data.append([x, y, z])
c1 = (
    Line3D()
    .add(
        "",
        data,
        xaxis3d_opts=opts.Axis3DOpts(Faker.clock, type_="value"),
        yaxis3d_opts=opts.Axis3DOpts(Faker.week_en, type_="value"),
        grid3d_opts=opts.Grid3DOpts(
            width=100, depth=100, rotate_speed=150, is_rotate=True
        ),
    )
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(
            max_=30, min_=0, range_color=Faker.visual_color
        ),
        title_opts=opts.TitleOpts(title="Line3D-旋转的弹簧"),
    )
)

# ---------------------------------------------- Line3d - Line3d_rectangular_projection ----------------------------------------------
import math

import pyecharts.options as opts
from pyecharts.charts import Line3D

week_en = "Saturday Friday Thursday Wednesday Tuesday Monday Sunday".split()
clock = (
    "12a 1a 2a 3a 4a 5a 6a 7a 8a 9a 10a 11a 12p "
    "1p 2p 3p 4p 5p 6p 7p 8p 9p 10p 11p".split()
)

data = []
for t in range(0, 25000):
    _t = t / 1000
    x = (1 + 0.25 * math.cos(75 * _t)) * math.cos(_t)
    y = (1 + 0.25 * math.cos(75 * _t)) * math.sin(_t)
    z = _t + 2.0 * math.sin(75 * _t)
    data.append([x, y, z])

c2 = (
    Line3D()
    .add(
        "",
        data,
        xaxis3d_opts=opts.Axis3DOpts(data=clock, type_="value"),
        yaxis3d_opts=opts.Axis3DOpts(data=week_en, type_="value"),
        grid3d_opts=opts.Grid3DOpts(width=100, height=100, depth=100),
    )
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(
            dimension=2,
            max_=30,
            min_=0,
            range_color=[
                "#313695",
                "#4575b4",
                "#74add1",
                "#abd9e9",
                "#e0f3f8",
                "#ffffbf",
                "#fee090",
                "#fdae61",
                "#f46d43",
                "#d73027",
                "#a50026",
            ],
        )
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="line3d",
    public=True,
)

swanlab.log(
    {
        "Line3d - Line3d_autorotate": c1,
        "Line3d - Line3d_rectangular_projection": c2
    }
)

