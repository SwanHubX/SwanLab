"""
// @author: ComPleHN
// @file: liquid.vue
// @time: 2025/5/27 18:45
// @description: 本文件是对于echarts的 水球图 图表测试 
"""
# ---------------------------------------------- Liquid - Liquid_shape_diamond ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid
from pyecharts.globals import SymbolType

c1 = (
    Liquid()
    .add("lq", [0.3, 0.7], is_outline_show=False, shape=SymbolType.DIAMOND)
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-Shape-Diamond"))
)

# ---------------------------------------------- Liquid - Multiple_liquid ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Grid, Liquid
from pyecharts.commons.utils import JsCode

l1 = (
    Liquid()
    .add("lq", [0.6, 0.7], center=["60%", "50%"])
    .set_global_opts(title_opts=opts.TitleOpts(title="多个 Liquid 显示"))
)

l2 = Liquid().add(
    "lq",
    [0.3254],
    center=["25%", "50%"],
    label_opts=opts.LabelOpts(
        font_size=50,
        formatter=JsCode(
            """function (param) {
                    return (Math.floor(param.value * 10000) / 100) + '%';
                }"""
        ),
        position="inside",
    ),
)

c2 = Grid().add(l1, grid_opts=opts.GridOpts()).add(l2, grid_opts=opts.GridOpts())

# ---------------------------------------------- Liquid - Liquid_shape_rect ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid
from pyecharts.globals import SymbolType

c3 = (
    Liquid()
    .add("lq", [0.3, 0.7], is_outline_show=False, shape=SymbolType.RECT)
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-Shape-rect"))
)

# ---------------------------------------------- Liquid - Liquid_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid

c4 = (
    Liquid()
    .add("lq", [0.6, 0.7])
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-基本示例"))
)

# ---------------------------------------------- Liquid - Liquid_data_precision ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid
from pyecharts.commons.utils import JsCode

c5 = (
    Liquid()
    .add(
        "lq",
        [0.3254],
        label_opts=opts.LabelOpts(
            font_size=50,
            formatter=JsCode(
                """function (param) {
                    return (Math.floor(param.value * 10000) / 100) + '%';
                }"""
            ),
            position="inside",
        ),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-数据精度"))
)

# ---------------------------------------------- Liquid - Liquid_without_outline ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid

c6 = (
    Liquid()
    .add("lq", [0.6, 0.7, 0.8], is_outline_show=False)
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-无边框"))
)

# ---------------------------------------------- Liquid - Liquid_shape_arrow ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Liquid
from pyecharts.globals import SymbolType

c7 = (
    Liquid()
    .add("lq", [0.3, 0.7], is_outline_show=False, shape=SymbolType.ARROW)
    .set_global_opts(title_opts=opts.TitleOpts(title="Liquid-Shape-arrow"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="liquid",
    public=True,
)

swanlab.log(
    {
        "Liquid - Liquid_shape_diamond": c1,
        "Liquid - Multiple_liquid": c2,
        "Liquid - Liquid_shape_rect": c3,
        "Liquid - Liquid_base": c4,
        "Liquid - Liquid_data_precision": c5,
        "Liquid - Liquid_without_outline": c6,
        "Liquid - Liquid_shape_arrow": c7
    }
)
