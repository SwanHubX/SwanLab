"""
// @author: ComPleHN
// @file: polar.vue
// @time: 2025/5/28 13:29
// @description: 本文件是对于echarts的 极坐标系 图表测试 
"""
# ---------------------------------------------- Polar - Polar_scatter_0 ----------------------------------------------
import random

from pyecharts import options as opts
from pyecharts.charts import Polar

data = [(i, random.randint(1, 100)) for i in range(101)]
c1 = (
    Polar()
    .add("", data, type_="scatter", label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-Scatter0"))
)

# ---------------------------------------------- Polar - Polar_love ----------------------------------------------
import math

from pyecharts import options as opts
from pyecharts.charts import Polar

data = []
for i in range(101):
    theta = i / 100 * 360
    r = 5 * (1 + math.sin(theta / 180 * math.pi))
    data.append([r, theta])
hour = [i for i in range(1, 25)]
c2 = (
    Polar()
    .add_schema(
        angleaxis_opts=opts.AngleAxisOpts(
            data=hour,
            type_="value",
            boundary_gap=False,
            start_angle=0,
            split_number=12,
            is_clockwise=True,
        )
    )
    .add("love", data, label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-Love"))
)

# ---------------------------------------------- Polar - Polar_radius ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Polar
from pyecharts.faker import Faker

c3 = (
    Polar()
    .add_schema(
        radiusaxis_opts=opts.RadiusAxisOpts(data=Faker.week, type_="category"),
        angleaxis_opts=opts.AngleAxisOpts(is_clockwise=True, max_=10),
    )
    .add("A", [1, 2, 3, 4, 3, 5, 1], type_="bar")
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-RadiusAxis"))
    .set_series_opts(label_opts=opts.LabelOpts(is_show=True))
)

# ---------------------------------------------- Polar - Two_value_axes_in_polar ----------------------------------------------
import math
import pyecharts.options as opts
from pyecharts.charts import Polar

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=line-polar

目前无法实现的功能:

1、赞无
"""

data = []

for i in range(0, 101):
    theta = i / 100 * 360
    r = 5 * (1 + math.sin(theta / 180 * math.pi))
    data.append([r, theta])

c4 = (
    Polar()
    .add(series_name="line", data=data, label_opts=opts.LabelOpts(is_show=False))
    .add_schema(
        angleaxis_opts=opts.AngleAxisOpts(
            start_angle=0, type_="value", is_clockwise=True
        )
    )
    .set_global_opts(
        tooltip_opts=opts.TooltipOpts(trigger="axis", axis_pointer_type="cross"),
        title_opts=opts.TitleOpts(title="极坐标双数值轴"),
    )
)

# ---------------------------------------------- Polar - Two_value_axes_in_polar_2 ----------------------------------------------
import math
import pyecharts.options as opts
from pyecharts.charts import Polar

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=line-polar2

目前无法实现的功能:

1、赞无
"""

data = []

for i in range(0, 360 + 1):
    t = i / 180 * math.pi
    r = math.sin(2 * t) * math.cos(2 * t)
    data.append([r, i])

c5 = (
    Polar()
    .add(
        series_name="line",
        data=data,
        label_opts=opts.LabelOpts(is_show=False),
        symbol_size=0,
    )
    .add_schema(
        angleaxis_opts=opts.AngleAxisOpts(
            start_angle=0, type_="value", is_clockwise=True
        ),
        radiusaxis_opts=opts.RadiusAxisOpts(min_=0),
    )
    .set_global_opts(
        tooltip_opts=opts.TooltipOpts(trigger="axis", axis_pointer_type="cross"),
        title_opts=opts.TitleOpts(title="极坐标双数值轴"),
    )
)

# ---------------------------------------------- Polar - Polar_scatter_1 ----------------------------------------------
import random

from pyecharts import options as opts
from pyecharts.charts import Polar

c6 = (
    Polar()
    .add("", [(10, random.randint(1, 100)) for i in range(300)], type_="scatter")
    .add("", [(11, random.randint(1, 100)) for i in range(300)], type_="scatter")
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-Scatter1"))
)

# ---------------------------------------------- Polar - Polar_angleaxis ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Polar
from pyecharts.faker import Faker

c7 = (
    Polar()
    .add_schema(angleaxis_opts=opts.AngleAxisOpts(data=Faker.week, type_="category"))
    .add("A", [1, 2, 3, 4, 3, 5, 1], type_="bar", stack="stack0")
    .add("B", [2, 4, 6, 1, 2, 3, 1], type_="bar", stack="stack0")
    .add("C", [1, 2, 3, 4, 1, 2, 5], type_="bar", stack="stack0")
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-AngleAxis"))
)

# ---------------------------------------------- Polar - Polar_flower ----------------------------------------------
import math

from pyecharts import options as opts
from pyecharts.charts import Polar

data = []
for i in range(361):
    t = i / 180 * math.pi
    r = math.sin(2 * t) * math.cos(2 * t)
    data.append([r, i])
c8 = (
    Polar()
    .add_schema(
        angleaxis_opts=opts.AngleAxisOpts(
            type_="value",
            boundary_gap=False,
            start_angle=0,
            split_number=12,
            is_clockwise=True,
        )
    )
    .add("flower", data, label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-Flower"))
)

# ---------------------------------------------- Polar - Polar_effectscatter ----------------------------------------------
import random

from pyecharts import options as opts
from pyecharts.charts import Polar


data = [(i, random.randint(1, 100)) for i in range(10)]
c9 = (
    Polar()
    .add(
        "",
        data,
        type_="effectScatter",
        effect_opts=opts.EffectOpts(scale=10, period=5),
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Polar-EffectScatter"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="polar",
    public=True,
)

swanlab.log(
    {
        "Polar - Polar_scatter_0": c1,
        "Polar - Polar_love": c2,
        "Polar - Polar_radius": c3,
        "Polar - Two_value_axes_in_polar": c4,
        "Polar - Two_value_axes_in_polar_2": c5,
        "Polar - Polar_scatter_1": c6,
        "Polar - Polar_angleaxis": c7,
        "Polar - Polar_flower": c8,
        "Polar - Polar_effectscatter": c9
    }
)
