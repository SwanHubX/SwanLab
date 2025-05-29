"""
// @author: ComPleHN
// @file: overlap.vue
// @time: 2025/5/28 12:45
// @description: 本文件是对于echarts的 层叠组件 图表测试 
"""
# ---------------------------------------------- Overlap - Overlap_bar_line ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar, Line
from pyecharts.faker import Faker

v1 = [2.0, 4.9, 7.0, 23.2, 25.6, 76.7, 135.6, 162.2, 32.6, 20.0, 6.4, 3.3]
v2 = [2.6, 5.9, 9.0, 26.4, 28.7, 70.7, 175.6, 182.2, 48.7, 18.8, 6.0, 2.3]
v3 = [2.0, 2.2, 3.3, 4.5, 6.3, 10.2, 20.3, 23.4, 23.0, 16.5, 12.0, 6.2]


c1= (
    Bar()
    .add_xaxis(Faker.months)
    .add_yaxis("蒸发量", v1)
    .add_yaxis("降水量", v2)
    .extend_axis(
        yaxis=opts.AxisOpts(
            axislabel_opts=opts.LabelOpts(formatter="{value} °C"), interval=5
        )
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Overlap-bar+line"),
        yaxis_opts=opts.AxisOpts(axislabel_opts=opts.LabelOpts(formatter="{value} ml")),
    )
)

line = Line().add_xaxis(Faker.months).add_yaxis("平均温度", v3, yaxis_index=1)
c1.overlap(line)

# ---------------------------------------------- Overlap - Overlap_line_scatter ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Line, Scatter
from pyecharts.faker import Faker

x = Faker.choose()
c2 = (
    Line()
    .add_xaxis(x)
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Overlap-line+scatter"))
)
scatter = (
    Scatter()
    .add_xaxis(x)
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
)
c2.overlap(scatter)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="overlap",
    public=True,
)

swanlab.log(
    {
        "Overlap - Overlap_bar_line": c1,
        "Overlap - Overlap_line_scatter": c2
    }
)
