"""
// @author: ComPleHN
// @file: radar.vue
// @time: 2025/5/29 10:19
// @description: 本文件是对于echarts的 雷达图 图表的测试
"""
# ---------------------------------------------- Radar - Basic_radar_chart ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Radar

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=radar

目前无法实现的功能:

1、雷达图周围的图例的 textStyle 暂时无法设置背景颜色
"""
v1 = [[4300, 10000, 28000, 35000, 50000, 19000]]
v2 = [[5000, 14000, 28000, 31000, 42000, 21000]]

c1 = (
    Radar(init_opts=opts.InitOpts(bg_color="#CCCCCC"))
    .add_schema(
        schema=[
            opts.RadarIndicatorItem(name="销售（sales）", max_=6500),
            opts.RadarIndicatorItem(name="管理（Administration）", max_=16000),
            opts.RadarIndicatorItem(name="信息技术（Information Technology）", max_=30000),
            opts.RadarIndicatorItem(name="客服（Customer Support）", max_=38000),
            opts.RadarIndicatorItem(name="研发（Development）", max_=52000),
            opts.RadarIndicatorItem(name="市场（Marketing）", max_=25000),
        ],
        splitarea_opt=opts.SplitAreaOpts(
            is_show=True, areastyle_opts=opts.AreaStyleOpts(opacity=1)
        ),
        textstyle_opts=opts.TextStyleOpts(color="#fff"),
    )
    .add(
        series_name="预算分配（Allocated Budget）",
        data=v1,
        linestyle_opts=opts.LineStyleOpts(color="#CD0000"),
    )
    .add(
        series_name="实际开销（Actual Spending）",
        data=v2,
        linestyle_opts=opts.LineStyleOpts(color="#5CACEE"),
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        title_opts=opts.TitleOpts(title="基础雷达图"), legend_opts=opts.LegendOpts()
    )
)

# ---------------------------------------------- Radar - Radar_selected_mode ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Radar

v1 = [[4300, 10000, 28000, 35000, 50000, 19000]]
v2 = [[5000, 14000, 28000, 31000, 42000, 21000]]
c2 = (
    Radar()
    .add_schema(
        schema=[
            opts.RadarIndicatorItem(name="销售", max_=6500),
            opts.RadarIndicatorItem(name="管理", max_=16000),
            opts.RadarIndicatorItem(name="信息技术", max_=30000),
            opts.RadarIndicatorItem(name="客服", max_=38000),
            opts.RadarIndicatorItem(name="研发", max_=52000),
            opts.RadarIndicatorItem(name="市场", max_=25000),
        ]
    )
    .add("预算分配", v1)
    .add("实际开销", v2)
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        legend_opts=opts.LegendOpts(selected_mode="single"),
        title_opts=opts.TitleOpts(title="Radar-单例模式"),
    )
)

# ---------------------------------------------- Radar - Radar_air_quality ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Radar

value_bj = [
    [55, 9, 56, 0.46, 18, 6, 1],
    [25, 11, 21, 0.65, 34, 9, 2],
    [56, 7, 63, 0.3, 14, 5, 3],
    [33, 7, 29, 0.33, 16, 6, 4],
    [42, 24, 44, 0.76, 40, 16, 5],
    [82, 58, 90, 1.77, 68, 33, 6],
    [74, 49, 77, 1.46, 48, 27, 7],
    [78, 55, 80, 1.29, 59, 29, 8],
    [267, 216, 280, 4.8, 108, 64, 9],
    [185, 127, 216, 2.52, 61, 27, 10],
    [39, 19, 38, 0.57, 31, 15, 11],
    [41, 11, 40, 0.43, 21, 7, 12],
]
value_sh = [
    [91, 45, 125, 0.82, 34, 23, 1],
    [65, 27, 78, 0.86, 45, 29, 2],
    [83, 60, 84, 1.09, 73, 27, 3],
    [109, 81, 121, 1.28, 68, 51, 4],
    [106, 77, 114, 1.07, 55, 51, 5],
    [109, 81, 121, 1.28, 68, 51, 6],
    [106, 77, 114, 1.07, 55, 51, 7],
    [89, 65, 78, 0.86, 51, 26, 8],
    [53, 33, 47, 0.64, 50, 17, 9],
    [80, 55, 80, 1.01, 75, 24, 10],
    [117, 81, 124, 1.03, 45, 24, 11],
    [99, 71, 142, 1.1, 62, 42, 12],
]
c_schema = [
    {"name": "AQI", "max": 300, "min": 5},
    {"name": "PM2.5", "max": 250, "min": 20},
    {"name": "PM10", "max": 300, "min": 5},
    {"name": "CO", "max": 5},
    {"name": "NO2", "max": 200},
    {"name": "SO2", "max": 100},
]
c3 = (
    Radar()
    .add_schema(schema=c_schema, shape="circle")
    .add("北京", value_bj, color="#f9713c")
    .add("上海", value_sh, color="#b3e4a1")
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Radar-空气质量"))
)

# ---------------------------------------------- Radar - Radar_angle_radius_axis ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Radar

data = [{"value": [4, -4, 2, 3, 0, 1], "name": "预算分配"}]
c_schema = [
    {"name": "销售", "max": 4, "min": -4},
    {"name": "管理", "max": 4, "min": -4},
    {"name": "技术", "max": 4, "min": -4},
    {"name": "客服", "max": 4, "min": -4},
    {"name": "研发", "max": 4, "min": -4},
    {"name": "市场", "max": 4, "min": -4},
]
c4 = (
    Radar()
    .set_colors(["#4587E7"])
    .add_schema(
        schema=c_schema,
        shape="circle",
        center=["50%", "50%"],
        radius="80%",
        angleaxis_opts=opts.AngleAxisOpts(
            min_=0,
            max_=360,
            is_clockwise=False,
            interval=5,
            axistick_opts=opts.AxisTickOpts(is_show=False),
            axislabel_opts=opts.LabelOpts(is_show=False),
            axisline_opts=opts.AxisLineOpts(is_show=False),
            splitline_opts=opts.SplitLineOpts(is_show=False),
        ),
        radiusaxis_opts=opts.RadiusAxisOpts(
            min_=-4,
            max_=4,
            interval=2,
            splitarea_opts=opts.SplitAreaOpts(
                is_show=True, areastyle_opts=opts.AreaStyleOpts(opacity=1)
            ),
        ),
        polar_opts=opts.PolarOpts(),
        splitarea_opt=opts.SplitAreaOpts(is_show=False),
        splitline_opt=opts.SplitLineOpts(is_show=False),
    )
    .add(
        series_name="预算",
        data=data,
        areastyle_opts=opts.AreaStyleOpts(opacity=0.1),
        linestyle_opts=opts.LineStyleOpts(width=1),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="radar",
    public=True,
)

swanlab.log(
    {
        "Radar - Basic_radar_chart": c1,
        "Radar - Radar_selected_mode": c2,
        "Radar - Radar_air_quality": c3,
        "Radar - Radar_angle_radius_axisr": c4
    }
)