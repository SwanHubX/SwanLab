"""
// @author: ComPleHN
// @file: effectscatter.vue
// @time: 2025/5/27 15:30
// @description: 本文件是对于echarts的 涟漪散点图 图表的测试
"""
# ---------------------------------------------- Effectscatter - Effectscatter_symbol ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import EffectScatter
from pyecharts.faker import Faker
from pyecharts.globals import SymbolType

c1 = (
    EffectScatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("", Faker.values(), symbol=SymbolType.ARROW)
    .set_global_opts(title_opts=opts.TitleOpts(title="EffectScatter-不同Symbol"))
)

# ---------------------------------------------- Effectscatter - Effectscatter_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import EffectScatter
from pyecharts.faker import Faker

c2 = (
    EffectScatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("", Faker.values())
    .set_global_opts(title_opts=opts.TitleOpts(title="EffectScatter-基本示例"))
)

# ---------------------------------------------- Effectscatter - Effectscatter_splitline ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import EffectScatter
from pyecharts.faker import Faker

c3 = (
    EffectScatter()
    .add_xaxis(Faker.choose())
    .add_yaxis("", Faker.values())
    .set_global_opts(
        title_opts=opts.TitleOpts(title="EffectScatter-显示分割线"),
        xaxis_opts=opts.AxisOpts(splitline_opts=opts.SplitLineOpts(is_show=True)),
        yaxis_opts=opts.AxisOpts(splitline_opts=opts.SplitLineOpts(is_show=True)),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="effect_scatter",
    public=True,
)

swanlab.log(
    {
        "Effectscatter - Effectscatter_symbol": c1,
        "Effectscatter - Effectscatter_base": c2,
        "Effectscatter - Effectscatter_splitline": c3
    }
)
