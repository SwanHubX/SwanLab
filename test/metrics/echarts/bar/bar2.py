"""
// @author: ComPleHN
// @file: bar2.vue
// @time: 2025/5/27 12:45
// @description: 由于bar类图表过多,我将分为多个文件进行测试
"""
# ---------------------------------------------- Bar - Bar_markline_type ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c1 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-MarkLine（指定类型）"))
    .set_series_opts(
        label_opts=opts.LabelOpts(is_show=False),
        markline_opts=opts.MarkLineOpts(
            data=[
                opts.MarkLineItem(type_="min", name="最小值"),
                opts.MarkLineItem(type_="max", name="最大值"),
                opts.MarkLineItem(type_="average", name="平均值"),
            ]
        ),
    )
)

# ---------------------------------------------- Bar - Bar_graphic_component ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.commons.utils import JsCode
from pyecharts.faker import Faker

c2 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Bar-Graphic 组件示例"),
        graphic_opts=[
            opts.GraphicGroup(
                graphic_item=opts.GraphicItem(
                    rotation=JsCode("Math.PI / 4"),
                    bounding="raw",
                    right=110,
                    bottom=110,
                    z=100,
                ),
                children=[
                    opts.GraphicRect(
                        graphic_item=opts.GraphicItem(
                            left="center", top="center", z=100
                        ),
                        graphic_shape_opts=opts.GraphicShapeOpts(width=400, height=50),
                        graphic_basicstyle_opts=opts.GraphicBasicStyleOpts(
                            fill="rgba(0,0,0,0.3)"
                        ),
                    ),
                    opts.GraphicText(
                        graphic_item=opts.GraphicItem(
                            left="center", top="center", z=100
                        ),
                        graphic_textstyle_opts=opts.GraphicTextStyleOpts(
                            text="pyecharts bar chart",
                            font="bold 26px Microsoft YaHei",
                            graphic_basicstyle_opts=opts.GraphicBasicStyleOpts(
                                fill="#fff"
                            ),
                        ),
                    ),
                ],
            )
        ],
    )
)

# ---------------------------------------------- Bar - Bar_waterfall_plot ----------------------------------------------
from pyecharts.charts import Bar
from pyecharts import options as opts

x_data = [f"11月{str(i)}日" for i in range(1, 12)]
y_total = [0, 900, 1245, 1530, 1376, 1376, 1511, 1689, 1856, 1495, 1292]
y_in = [900, 345, 393, "-", "-", 135, 178, 286, "-", "-", "-"]
y_out = ["-", "-", "-", 108, 154, "-", "-", "-", 119, 361, 203]


c3 = (
    Bar()
    .add_xaxis(xaxis_data=x_data)
    .add_yaxis(
        series_name="",
        y_axis=y_total,
        stack="总量",
        itemstyle_opts=opts.ItemStyleOpts(color="rgba(0,0,0,0)"),
    )
    .add_yaxis(series_name="收入", y_axis=y_in, stack="总量")
    .add_yaxis(series_name="支出", y_axis=y_out, stack="总量")
    .set_global_opts(yaxis_opts=opts.AxisOpts(type_="value"))
)

# ---------------------------------------------- Bar - Bar_border_radius ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.commons.utils import JsCode
from pyecharts.faker import Faker

c4 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values(), category_gap="60%")
    .set_series_opts(
        itemstyle_opts={
            "normal": {
                "color": JsCode(
                    """new echarts.graphic.LinearGradient(0, 0, 0, 1, [{
                offset: 0,
                color: 'rgba(0, 244, 255, 1)'
            }, {
                offset: 1,
                color: 'rgba(0, 77, 167, 1)'
            }], false)"""
                ),
                "barBorderRadius": [30, 30, 30, 30],
                "shadowColor": "rgb(0, 160, 221)",
            }
        }
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-渐变圆柱"))
)

# ---------------------------------------------- Bar - Bar_same_series_gap ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker


c5 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values(), category_gap="80%")
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-单系列柱间距离"))
)

# ---------------------------------------------- Bar - Bar_datazoom_inside ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c6 = (
    Bar()
    .add_xaxis(Faker.days_attrs)
    .add_yaxis("商家A", Faker.days_values, color=Faker.rand_color())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Bar-DataZoom（inside）"),
        datazoom_opts=opts.DataZoomOpts(type_="inside"),
    )
)

# ---------------------------------------------- Bar - Bar_is_selected ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker


c7 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Bar-默认取消显示某 Series"),
        legend_opts=opts.LegendOpts(selected_map={"商家B": False}),
    )
)

# ---------------------------------------------- Bar - Bar_reversal_axis ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c8 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .reversal_axis()
    .set_series_opts(label_opts=opts.LabelOpts(position="right"))
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-翻转 XY 轴"))
)

# ---------------------------------------------- Bar - Bar_markpoint_custom ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

x, y = Faker.choose(), Faker.values()
c9 = (
    Bar()
    .add_xaxis(x)
    .add_yaxis(
        "商家A",
        y,
        markpoint_opts=opts.MarkPointOpts(
            data=[opts.MarkPointItem(name="自定义标记点", coord=[x[2], y[2]], value=y[2])]
        ),
    )
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-MarkPoint（自定义）"))
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
)

# ---------------------------------------------- Bar - Bar_base_with_animation ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker


c10 = (
    Bar(
        init_opts=opts.InitOpts(
            animation_opts=opts.AnimationOpts(
                animation_delay=1000, animation_easing="elasticOut"
            )
        )
    )
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-动画配置基本示例", subtitle="我是副标题"))
)

# ---------------------------------------------- Bar - Bar_histogram ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c11 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values(), category_gap=0, color=Faker.rand_color())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-直方图"))
)

# ---------------------------------------------- Bar - Bar_markline_custom ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c12 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-MarkLine（自定义）"))
    .set_series_opts(
        label_opts=opts.LabelOpts(is_show=False),
        markline_opts=opts.MarkLineOpts(
            data=[opts.MarkLineItem(y=50, name="yAxis=50")]
        ),
    )
)

# ---------------------------------------------- Bar - Bar_base ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker


c13 = (
    Bar()
    .add_xaxis(Faker.choose())
    .add_yaxis("商家A", Faker.values())
    .add_yaxis("商家B", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="Bar-基本示例", subtitle="我是副标题"))
)

# ---------------------------------------------- Bar - Bar_datazoom_both ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar
from pyecharts.faker import Faker

c14 = (
    Bar()
    .add_xaxis(Faker.days_attrs)
    .add_yaxis("商家A", Faker.days_values, color=Faker.rand_color())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Bar-DataZoom（slider+inside）"),
        datazoom_opts=[opts.DataZoomOpts(), opts.DataZoomOpts(type_="inside")],
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="bar2",
    public=True,
)

swanlab.log(
    {
        "Bar_markline_type": c1,
        "Bar_graphic_component": c2,
        "Bar_waterfall_plot": c3,
        "Bar_border_radius": c4,
        "Bar_same_series_gap": c5,
        "Bar_datazoom_inside": c6,
        "Bar_is_selected": c7,
        "Bar_reversal_axis": c8,
        "Bar_markpoint_custom": c9,
        "Bar_base_with_animation": c10,
        "Bar_histogram": c11,
        "Bar_markline_custom": c12,
        "Bar_base": c13,
        "Bar_datazoom_both": c14
        
    }
)