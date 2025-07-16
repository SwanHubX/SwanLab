"""
// @author: ComPleHN
// @file: scatter.vue
// @time: 2025/5/29 10:32
// @description: 本文件是对于echarts的 散点图 图表的测试
"""
# ---------------------------------------------- Scatter - Basic_scatter_chart ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Scatter

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=scatter-simple

目前无法实现的功能:

1、暂无
"""

data = [
    [10.0, 8.04],
    [8.0, 6.95],
    [13.0, 7.58],
    [9.0, 8.81],
    [11.0, 8.33],
    [14.0, 9.96],
    [6.0, 7.24],
    [4.0, 4.26],
    [12.0, 10.84],
    [7.0, 4.82],
    [5.0, 5.68],
]
data.sort(key=lambda x: x[0])
x_data = [d[0] for d in data]
y_data = [d[1] for d in data]

c1 = (
    Scatter()
    .add_xaxis(xaxis_data=x_data)
    .add_yaxis(
        series_name="",
        y_axis=y_data,
        symbol_size=20,
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_series_opts()
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value", splitline_opts=opts.SplitLineOpts(is_show=True)
        ),
        yaxis_opts=opts.AxisOpts(
            type_="value",
            axistick_opts=opts.AxisTickOpts(is_show=True),
            splitline_opts=opts.SplitLineOpts(is_show=True),
        ),
        tooltip_opts=opts.TooltipOpts(is_show=False),
    )
)

# ---------------------------------------------- Scatter - Scatter_multi_dimension ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Scatter
from pyecharts.commons.utils import JsCode
from pyecharts.faker import Faker

c2 = (
    Scatter()
    .add_xaxis(Faker.choose())
    .add_yaxis(
        "商家A",
        [list(z) for z in zip(Faker.values(), Faker.choose())],
        label_opts=opts.LabelOpts(
            formatter=JsCode(
                "function(params){return params.value[1] +' : '+ params.value[2];}"
            )
        ),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Scatter-多维度数据"),
        tooltip_opts=opts.TooltipOpts(
            formatter=JsCode(
                "function (params) {return params.name + ' : ' + params.value[2];}"
            )
        ),
        visualmap_opts=opts.VisualMapOpts(
            type_="color", max_=150, min_=20, dimension=1
        ),
    )
)

# ----------------------------------------------Scatter - Scatter_visualmap_color  ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Scatter
from pyecharts.faker import Faker

c3 = (
    Scatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Scatter-VisualMap(Color)"),
        visualmap_opts=opts.VisualMapOpts(max_=150),
    )
)

# ---------------------------------------------- Scatter - Scatter_splitline ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Scatter
from pyecharts.faker import Faker

c4 = (
    Scatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Scatter-显示分割线"),
        xaxis_opts=opts.AxisOpts(splitline_opts=opts.SplitLineOpts(is_show=True)),
        yaxis_opts=opts.AxisOpts(splitline_opts=opts.SplitLineOpts(is_show=True)),
    )
)

# ---------------------------------------------- Scatter - Scatter_visualmap_size ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Scatter
from pyecharts.faker import Faker

c5 = (
    Scatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Scatter-VisualMap(Size)"),
        visualmap_opts=opts.VisualMapOpts(type_="size", max_=150, min_=20),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="scatter",
    public=True,
)

swanlab.log(
    {
        "Scatter - Basic_scatter_chart": c1,
        "Scatter - Scatter_multi_dimension": c2,
        "Scatter - Scatter_visualmap_color": c3,
        "Scatter - Scatter_splitline": c4,
        "Scatter - Scatter_visualmap_size": c5
    }
)
