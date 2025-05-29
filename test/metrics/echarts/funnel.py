"""
// @author: ComPleHN
// @file: funnel.vue
// @time: 2025/5/27 15:36
// @description: 本文件是对于echarts的 漏斗图 图表测试
"""
# ---------------------------------------------- Funnel - Funnel_sort_ascending ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Funnel
from pyecharts.faker import Faker

c1 = (
    Funnel()
    .add(
        "商品",
        [list(z) for z in zip(Faker.choose(), Faker.values())],
        sort_="ascending",
        label_opts=opts.LabelOpts(position="inside"),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Funnel-Sort（ascending）"))
)

# ---------------------------------------------- Funnel - Funnel_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Funnel
from pyecharts.faker import Faker

c2 = (
    Funnel()
    .add("商品", [list(z) for z in zip(Faker.choose(), Faker.values())])
    .set_global_opts(title_opts=opts.TitleOpts(title="Funnel-基本示例"))
)

# ---------------------------------------------- Funnel - Funnel_chart ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import Funnel

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=funnel

目前无法实现的功能:

1、暂时无法对漏斗图的长宽等范围操作进行修改
"""
x_data = ["展现", "点击", "访问", "咨询", "订单"]
y_data = [100, 80, 60, 40, 20]

data = [[x_data[i], y_data[i]] for i in range(len(x_data))]

c3 = (
    Funnel()
    .add(
        series_name="",
        data_pair=data,
        gap=2,
        tooltip_opts=opts.TooltipOpts(trigger="item", formatter="{a} <br/>{b} : {c}%"),
        label_opts=opts.LabelOpts(is_show=True, position="inside"),
        itemstyle_opts=opts.ItemStyleOpts(border_color="#fff", border_width=1),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="漏斗图", subtitle="纯属虚构"))
)

# ---------------------------------------------- Funnel - Funnel_label_inside ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Funnel
from pyecharts.faker import Faker


c4 = (
    Funnel()
    .add(
        "商品",
        [list(z) for z in zip(Faker.choose(), Faker.values())],
        label_opts=opts.LabelOpts(position="inside"),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Funnel-Label（inside)"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="funnel",
    public=True,
)

swanlab.log(
    {
        "Funnel - Funnel_sort_ascending": c1,
        "Funnel - Funnel_base": c2,
        "Funnel - Funnel_chart": c3,
        "Funnel - Funnel_label_inside": c4
    }
)
