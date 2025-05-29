"""
// @author: ComPleHN
// @file: mapglobe.vue
// @time: 2025/5/28 12:41
// @description: 本文件是对于echarts的 Globe地图 图表测试 
"""
# ---------------------------------------------- Mapglobe - Map_globe_base ----------------------------------------------
import pyecharts.options as opts
from pyecharts.charts import MapGlobe
from pyecharts.faker import POPULATION

data = [x for _, x in POPULATION[1:]]
low, high = min(data), max(data)

c1 = (
    MapGlobe()
    .add_schema()
    .add(
        maptype="world",
        series_name="World Population",
        data_pair=POPULATION[1:],
        is_map_symbol_show=False,
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        visualmap_opts=opts.VisualMapOpts(
            min_=low,
            max_=high,
            range_text=["max", "min"],
            is_calculable=True,
            range_color=["lightskyblue", "yellow", "orangered"],
        )
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="liquid",
    public=True,
)

swanlab.log(
    {
        "Mapglobe - Map_globe_base": c1
    }
)
