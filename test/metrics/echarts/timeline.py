"""
// @author: ComPleHN
// @file: timeline.vue
// @time: 2025/5/29 10:52
// @description: 本文件是对于echarts的 时间轴组件 图表的测试
"""
# ---------------------------------------------- Timeline - Timeline_pie ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Pie, Timeline
from pyecharts.faker import Faker

attr = Faker.choose()
c1 = Timeline()
for i in range(2015, 2020):
    pie = (
        Pie()
        .add(
            "商家A",
            [list(z) for z in zip(attr, Faker.values())],
            rosetype="radius",
            radius=["30%", "55%"],
        )
        .set_global_opts(title_opts=opts.TitleOpts("某商店{}年营业额".format(i)))
    )
    c1.add(pie, "{}年".format(i))

# ---------------------------------------------- Timeline - Timeline_bar_with_graphic ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar, Timeline
from pyecharts.commons.utils import JsCode
from pyecharts.faker import Faker

x = Faker.choose()
c2 = Timeline()
for i in range(2015, 2020):
    bar = (
        Bar()
        .add_xaxis(x)
        .add_yaxis("商家A", Faker.values())
        .add_yaxis("商家B", Faker.values())
        .set_global_opts(
            title_opts=opts.TitleOpts("某商店{}年营业额 - With Graphic 组件".format(i)),
            graphic_opts=[
                opts.GraphicGroup(
                    graphic_item=opts.GraphicItem(
                        rotation=JsCode("Math.PI / 4"),
                        bounding="raw",
                        right=100,
                        bottom=110,
                        z=100,
                    ),
                    children=[
                        opts.GraphicRect(
                            graphic_item=opts.GraphicItem(
                                left="center", top="center", z=100
                            ),
                            graphic_shape_opts=opts.GraphicShapeOpts(
                                width=400, height=50
                            ),
                            graphic_basicstyle_opts=opts.GraphicBasicStyleOpts(
                                fill="rgba(0,0,0,0.3)"
                            ),
                        ),
                        opts.GraphicText(
                            graphic_item=opts.GraphicItem(
                                left="center", top="center", z=100
                            ),
                            graphic_textstyle_opts=opts.GraphicTextStyleOpts(
                                text="某商店{}年营业额".format(i),
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
    c2.add(bar, "{}年".format(i))

# ---------------------------------------------- Timeline - Timeline_bmap ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import BMap, Timeline
from pyecharts.faker import Faker

c3 = Timeline()
c3.add_schema(pos_left="50%", pos_right="10px", pos_bottom="15px")
for i in range(2015, 2020):
    bmap = (
        BMap()
        .add_schema(baidu_ak="FAKE_AK", center=[120.13066322374, 30.240018034923])
        .add(
            "bmap",
            [list(z) for z in zip(Faker.provinces, Faker.values())],
            type_="heatmap",
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="Timeline-BMap-热力图-{}年".format(i)),
            visualmap_opts=opts.VisualMapOpts(pos_bottom="center", pos_right="10px"),
            tooltip_opts=opts.TooltipOpts(formatter=None),
        )
    )
    c3.add(bmap, "{}年".format(i))

# ---------------------------------------------- Timeline - Timeline_sankey ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Sankey, Timeline
from pyecharts.faker import Faker

c4 = Timeline()
names = ("商家A", "商家B", "商家C")
nodes = [{"name": name} for name in names]
for i in range(2015, 2020):
    links = [
        {"source": names[0], "target": names[1], "value": Faker.values()[0]},
        {"source": names[1], "target": names[2], "value": Faker.values()[0]},
    ]
    sankey = (
        Sankey()
        .add(
            "sankey",
            nodes,
            links,
            linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),
            label_opts=opts.LabelOpts(position="right"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="{}年商店（A, B, C）营业额差".format(i))
        )
    )
    c4.add(sankey, "{}年".format(i))
    
# ---------------------------------------------- Timeline - Timeline_bar_reversal ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar, Timeline
from pyecharts.faker import Faker

c5 = Timeline()
for i in range(2015, 2020):
    bar = (
        Bar()
        .add_xaxis(Faker.choose())
        .add_yaxis("商家A", Faker.values(), label_opts=opts.LabelOpts(position="right"))
        .add_yaxis("商家B", Faker.values(), label_opts=opts.LabelOpts(position="right"))
        .reversal_axis()
        .set_global_opts(
            title_opts=opts.TitleOpts("Timeline-Bar-Reversal (时间: {} 年)".format(i))
        )
    )
    c5.add(bar, "{}年".format(i))

# ---------------------------------------------- Timeline - Timeline_bar ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar, Timeline
from pyecharts.faker import Faker

x = Faker.choose()
c6 = Timeline()
for i in range(2015, 2020):
    bar = (
        Bar()
        .add_xaxis(x)
        .add_yaxis("商家A", Faker.values())
        .add_yaxis("商家B", Faker.values())
        .set_global_opts(title_opts=opts.TitleOpts("某商店{}年营业额".format(i)))
    )
    c6.add(bar, "{}年".format(i))

# ---------------------------------------------- Timeline - Timeline_multi_axis ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar, Timeline
from pyecharts.faker import Faker

c7 = Timeline()
for i in range(2015, 2020):
    bar = (
        Bar()
        .add_xaxis(Faker.choose())
        .add_yaxis("商家A", Faker.values())
        .add_yaxis("商家B", Faker.values())
        .set_global_opts(title_opts=opts.TitleOpts("某商店{}年营业额".format(i)))
    )
    c7.add(bar, "{}年".format(i))
    
# ---------------------------------------------- Timeline - Timeline_map ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Map, Timeline
from pyecharts.faker import Faker

c8 = Timeline()
for i in range(2015, 2020):
    map0 = (
        Map()
        .add("商家A", [list(z) for z in zip(Faker.provinces, Faker.values())], "china")
        .set_global_opts(
            title_opts=opts.TitleOpts(title="Map-{}年某些数据".format(i)),
            visualmap_opts=opts.VisualMapOpts(max_=200),
        )
    )
    c8.add(map0, "{}年".format(i))

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="timeline",
    public=True,
)

swanlab.log(
    {
        "Timeline - Timeline_pie": c1,
        "Timeline - Timeline_bar_with_graphic": c2,
        "Timeline - Timeline_bmap": c3,
        "Timeline - Timeline_sankey": c4,
        "Timeline - Timeline_bar_reversal": c5,
        "Timeline - Timeline_bar": c6,
        "Timeline - Timeline_multi_axis": c7,
        "Timeline - Timeline_map": c8
    }
)
