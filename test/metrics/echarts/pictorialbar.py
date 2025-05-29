"""
// @author: ComPleHN
// @file: pictorialbar.vue
// @time: 2025/5/28 12:58
// @description: 本文件是对于echarts的 象形柱图 图表测试,本文件有注意事项请看本目录下 README.md
"""
# ---------------------------------------------- Pictorialbar - Pictorialbar_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import PictorialBar
from pyecharts.globals import SymbolType

location = ["山西", "四川", "西藏", "北京", "上海", "内蒙古", "云南", "黑龙江", "广东", "福建"]
values = [13, 42, 67, 81, 86, 94, 166, 220, 249, 262]

c1 = (
    PictorialBar()
    .add_xaxis(location)
    .add_yaxis(
        "",
        values,
        label_opts=opts.LabelOpts(is_show=False),
        symbol_size=18,
        symbol_repeat="fixed",
        symbol_offset=[0, 0],
        is_symbol_clip=True,
        symbol=SymbolType.ROUND_RECT,
    )
    .reversal_axis()
    .set_global_opts(
        title_opts=opts.TitleOpts(title="PictorialBar-各省份人口数量（虚假数据）"),
        xaxis_opts=opts.AxisOpts(is_show=False),
        yaxis_opts=opts.AxisOpts(
            axistick_opts=opts.AxisTickOpts(is_show=False),
            axisline_opts=opts.AxisLineOpts(
                linestyle_opts=opts.LineStyleOpts(opacity=0)
            ),
        ),
    )
)

# ---------------------------------------------- Pictorialbar - Pictorialbar_multi_custom_symbols ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import PictorialBar

location = ["山西", "四川", "西藏", "北京", "上海", "内蒙古", "云南", "黑龙江", "广东", "福建"]
values = [13, 42, 67, 81, 86, 94, 166, 220, 249, 262]


with open("../assets/echarts/pictorialbar/symbol.json", "r", encoding="utf-8") as f:
    symbols = json.load(f)

c2 = (
    PictorialBar()
    .add_xaxis(["reindeer", "ship", "plane", "train", "car"])
    .add_yaxis(
        "2015",
        [
            {"value": 157, "symbol": symbols["reindeer"]},
            {"value": 21, "symbol": symbols["ship"]},
            {"value": 66, "symbol": symbols["plane"]},
            {"value": 78, "symbol": symbols["train"]},
            {"value": 123, "symbol": symbols["car"]},
        ],
        label_opts=opts.LabelOpts(is_show=False),
        symbol_size=22,
        symbol_repeat="fixed",
        symbol_offset=[0, 5],
        is_symbol_clip=True,
    )
    .add_yaxis(
        "2016",
        [
            {"value": 184, "symbol": symbols["reindeer"]},
            {"value": 29, "symbol": symbols["ship"]},
            {"value": 73, "symbol": symbols["plane"]},
            {"value": 91, "symbol": symbols["train"]},
            {"value": 95, "symbol": symbols["car"]},
        ],
        label_opts=opts.LabelOpts(is_show=False),
        symbol_size=22,
        symbol_repeat="fixed",
        symbol_offset=[0, -25],
        is_symbol_clip=True,
    )
    .reversal_axis()
    .set_global_opts(
        title_opts=opts.TitleOpts(title="PictorialBar-Vehicles in X City"),
        xaxis_opts=opts.AxisOpts(is_show=False),
        yaxis_opts=opts.AxisOpts(
            axistick_opts=opts.AxisTickOpts(is_show=False),
            axisline_opts=opts.AxisLineOpts(
                linestyle_opts=opts.LineStyleOpts(opacity=0)
            ),
        ),
    )
)

# ---------------------------------------------- Pictorialbar - Pictorialbar_custom_symbol ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import PictorialBar

location = ["山西", "四川", "西藏", "北京", "上海", "内蒙古", "云南", "黑龙江", "广东", "福建"]
values = [13, 42, 67, 81, 86, 94, 166, 220, 249, 262]


with open("../assets/echarts/pictorialbar/symbol.json", "r", encoding="utf-8") as f:
    symbols = json.load(f)

c3 = (
    PictorialBar()
    .add_xaxis(location)
    .add_yaxis(
        "",
        values,
        label_opts=opts.LabelOpts(is_show=False),
        symbol_size=22,
        symbol_repeat="fixed",
        symbol_offset=[0, -5],
        is_symbol_clip=True,
        symbol=symbols["boy"],
    )
    .reversal_axis()
    .set_global_opts(
        title_opts=opts.TitleOpts(title="PictorialBar-自定义 Symbol"),
        xaxis_opts=opts.AxisOpts(is_show=False),
        yaxis_opts=opts.AxisOpts(
            axistick_opts=opts.AxisTickOpts(is_show=False),
            axisline_opts=opts.AxisLineOpts(
                linestyle_opts=opts.LineStyleOpts(opacity=0)
            ),
        ),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="pictorialbar",
    public=True,
)

swanlab.log(
    {
        "Pictorialbar - Pictorialbar_base": c1,
        "Pictorialbar - Pictorialbar_multi_custom_symbols": c2,
        "Pictorialbar - Pictorialbar_custom_symbol": c3
    }
)
