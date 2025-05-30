"""
// @author: ComPleHN
// @file: gauge.vue
// @time: 2025/5/27 15:47
// @description: 本文件是对于echarts的 仪表盘 图表测试
"""
# ---------------------------------------------- Gauge - Gauge ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Gauge

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=gauge

目前无法实现的功能:

1、暂无
"""

c1 = (
    Gauge()
    .add(series_name="业务指标", data_pair=[["完成率", 55.5]])
    .set_global_opts(
        legend_opts=opts.LegendOpts(is_show=False),
        tooltip_opts=opts.TooltipOpts(is_show=True, formatter="{a} <br/>{b} : {c}%"),
    )
)

# ---------------------------------------------- Gauge - Gauge_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Gauge

c2 = (
    Gauge()
    .add("", [("完成率", 66.6)])
    .set_global_opts(title_opts=opts.TitleOpts(title="Gauge-基本示例"))
)

# ---------------------------------------------- Gauge - Gauge_color ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Gauge

c3 = (
    Gauge()
    .add(
        "业务指标",
        [("完成率", 55.5)],
        axisline_opts=opts.AxisLineOpts(
            linestyle_opts=opts.LineStyleOpts(
                color=[(0.3, "#67e0e3"), (0.7, "#37a2da"), (1, "#fd666d")], width=30
            )
        ),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Gauge-不同颜色"),
        legend_opts=opts.LegendOpts(is_show=False),
    )
)

# ---------------------------------------------- Gauge - Gauge_label_title_setting ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Gauge

c4 = (
    Gauge()
    .add(
        "",
        [("完成率", 66.6)],
        title_label_opts=opts.LabelOpts(
            font_size=40, color="blue", font_family="Microsoft YaHei"
        ),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Gauge-改变轮盘内的字体"))
)

# ---------------------------------------------- Gauge - Gauge_change_color ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Gauge

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://gallery.echartsjs.com/editor.html?c=xH1vxib94f

目前无法实现的功能:

1、暂无
"""

c5 = (
    Gauge()
    .add(series_name="业务指标", data_pair=[["完成率", 55.5]])
    .set_global_opts(
        legend_opts=opts.LegendOpts(is_show=False),
        tooltip_opts=opts.TooltipOpts(is_show=True, formatter="{a} <br/>{b} : {c}%"),
    )
    .set_series_opts(
        axisline_opts=opts.AxisLineOpts(
            linestyle_opts=opts.LineStyleOpts(
                color=[[0.3, "#67e0e3"], [0.7, "#37a2da"], [1, "#fd666d"]], width=30
            )
        )
    )
)

# ---------------------------------------------- Gauge - Gauge_splitnum_label ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Gauge

c6 = (
    Gauge()
    .add(
        "业务指标",
        [("完成率", 55.5)],
        split_number=5,
        axisline_opts=opts.AxisLineOpts(
            linestyle_opts=opts.LineStyleOpts(
                color=[(0.3, "#67e0e3"), (0.7, "#37a2da"), (1, "#fd666d")], width=30
            )
        ),
        detail_label_opts=opts.LabelOpts(formatter="{value}"),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Gauge-分割段数-Label"),
        legend_opts=opts.LegendOpts(is_show=False),
    )
)

# ---------------------------------------------- Gauge - Gauge_change_radius ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Gauge

c7 = (
    Gauge()
    .add("", [("完成率", 66.6)], radius="50%")
    .set_global_opts(title_opts=opts.TitleOpts(title="Gauge-修改 Radius 为 50%"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="gauge",
    public=True,
)

swanlab.log(
    {
        "Gauge - Gauge": c1,
        "Gauge - Gauge_base": c2,
        "Gauge - Gauge_color": c3,
        "Gauge - Gauge_label_title_setting": c4,
        "Gauge - Gauge_change_color": c5,
        "Gauge - Gauge_splitnum_label": c6,
        "Gauge - Gauge_change_radius": c7

    }
)