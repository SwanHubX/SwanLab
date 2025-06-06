"""
// @author: ComPleHN
// @file: feo.vue
// @time: 2025/5/27 15:53
// @description: 本文件是对于echarts的 地理坐标 图表测试
"""
# ---------------------------------------------- Geo - Geo_chart_countries_js ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.datasets import register_url

# try:
#     register_url("https://assets.pyecharts.org/assets/maps/")
# except Exception:
#     import ssl
#
#     ssl._create_default_https_context = ssl._create_unverified_context
#     register_url("https://assets.pyecharts.org/assets/maps/")

c1 = (
    Geo()
    .add_schema(maptype="瑞士")
    .set_global_opts(title_opts=opts.TitleOpts(title="瑞士"))
)

# ---------------------------------------------- Geo - Geo_lines_background ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.globals import ChartType, SymbolType

c2 = (
    Geo()
    .add_schema(
        maptype="china",
        itemstyle_opts=opts.ItemStyleOpts(color="#323c48", border_color="#111"),
    )
    .add(
        "",
        [("广州", 55), ("北京", 66), ("杭州", 77), ("重庆", 88)],
        type_=ChartType.EFFECT_SCATTER,
        color="white",
    )
    .add(
        "geo",
        [("广州", "上海"), ("广州", "北京"), ("广州", "杭州"), ("广州", "重庆")],
        type_=ChartType.LINES,
        effect_opts=opts.EffectOpts(
            symbol=SymbolType.ARROW, symbol_size=6, color="blue"
        ),
        linestyle_opts=opts.LineStyleOpts(curve=0.2),
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Geo-Lines-background"))
)

# ---------------------------------------------- Geo - Geo_visualmap_piecewise ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.faker import Faker

c3 = (
    Geo()
    .add_schema(maptype="china")
    .add("geo", [list(z) for z in zip(Faker.provinces, Faker.values())])
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(is_piecewise=True),
        title_opts=opts.TitleOpts(title="Geo-VisualMap（分段型）"),
    )
)

# ---------------------------------------------- Geo - Geo_lines ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.globals import ChartType, SymbolType

c4 = (
    Geo()
    .add_schema(maptype="china")
    .add(
        "",
        [("广州", 55), ("北京", 66), ("杭州", 77), ("重庆", 88)],
        type_=ChartType.EFFECT_SCATTER,
        color="white",
    )
    .add(
        "geo",
        [("广州", "上海"), ("广州", "北京"), ("广州", "杭州"), ("广州", "重庆")],
        type_=ChartType.LINES,
        effect_opts=opts.EffectOpts(
            symbol=SymbolType.ARROW, symbol_size=6, color="blue"
        ),
        linestyle_opts=opts.LineStyleOpts(curve=0.2),
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Geo-Lines"))
)

# ---------------------------------------------- Geo - Geo_guangdong ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.faker import Faker
from pyecharts.globals import ChartType

c5 = (
    Geo()
    .add_schema(maptype="广东")
    .add(
        "geo",
        [list(z) for z in zip(Faker.guangdong_city, Faker.values())],
        type_=ChartType.HEATMAP,
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(), title_opts=opts.TitleOpts(title="Geo-广东地图")
    )
)

# ---------------------------------------------- Geo - Geo_heatmap ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.faker import Faker
from pyecharts.globals import ChartType

c6 = (
    Geo()
    .add_schema(maptype="china")
    .add(
        "geo",
        [list(z) for z in zip(Faker.provinces, Faker.values())],
        type_=ChartType.HEATMAP,
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(),
        title_opts=opts.TitleOpts(title="Geo-HeatMap"),
    )
)

# ---------------------------------------------- Geo - Geo_effectscatter ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.faker import Faker
from pyecharts.globals import ChartType

c7 = (
    Geo()
    .add_schema(maptype="china")
    .add(
        "geo",
        [list(z) for z in zip(Faker.provinces, Faker.values())],
        type_=ChartType.EFFECT_SCATTER,
    )
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(title_opts=opts.TitleOpts(title="Geo-EffectScatter"))
)

# ---------------------------------------------- Geo - Geo_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.faker import Faker

c8 = (
    Geo()
    .add_schema(maptype="china")
    .add("geo", [list(z) for z in zip(Faker.provinces, Faker.values())])
    .set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(), title_opts=opts.TitleOpts(title="Geo-基本示例")
    )
)

# ---------------------------------------------- Geo - Geo_echart_china_js ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Geo
from pyecharts.datasets import register_url

# try:
#     register_url("https://assets.pyecharts.org/assets/v5/maps/")
# except Exception:
#     import ssl
#
#     ssl._create_default_https_context = ssl._create_unverified_context
#     register_url("https://assets.pyecharts.org/assets/v5/maps/")

c9 = (
    Geo()
    .add_schema(maptype="海淀区")
    .set_global_opts(title_opts=opts.TitleOpts(title="海淀区"))
)




import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="geo",
    public=True,
)

swanlab.log(
    {
        "Geo - Geo_chart_countries_js": c1,
        "Geo - Geo_lines_background": c2,
        "Geo - Geo_visualmap_piecewise": c3,
        "Geo - Geo_lines": c4,
        "Geo - Geo_guangdong": c5,
        "Geo - Geo_heatmap": c6,
        "Geo - Geo_effectscatter": c7,
        "Geo - Geo_base": c8,
        "Geo - Geo_echart_china_js": c9
    }
)
