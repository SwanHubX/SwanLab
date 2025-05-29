"""
// @author: ComPleHN
// @file: sankey.vue
// @time: 2025/5/29 10:24
// @description: 本文件是对于echarts的 桑葚图 图表的测试,本文件有注意事项请看本目录下 README.md
"""
# ---------------------------------------------- Sankey - Sankey_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Sankey

nodes = [
    {"name": "category1"},
    {"name": "category2"},
    {"name": "category3"},
    {"name": "category4"},
    {"name": "category5"},
    {"name": "category6"},
]

links = [
    {"source": "category1", "target": "category2", "value": 10},
    {"source": "category2", "target": "category3", "value": 15},
    {"source": "category3", "target": "category4", "value": 20},
    {"source": "category5", "target": "category6", "value": 25},
]
c1 = (
    Sankey()
    .add(
        "sankey",
        nodes,
        links,
        linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),
        label_opts=opts.LabelOpts(position="right"),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Sankey-基本示例"))
)

# ---------------------------------------------- Sankey - Sankey_diagram ----------------------------------------------
import asyncio
from aiohttp import TCPConnector, ClientSession

import pyecharts.options as opts
from pyecharts.charts import Sankey

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=sankey-energy

目前无法实现的功能:

1、label 和图的层次有点问题
"""


async def get_json_data(url: str) -> dict:
    async with ClientSession(connector=TCPConnector(ssl=False)) as session:
        async with session.get(url=url) as response:
            return await response.json()


# 获取官方的数据
data = asyncio.run(
    get_json_data(url="https://echarts.apache.org/examples/data/asset/data/energy.json")
)

c2 = (
    Sankey()
    .add(
        series_name="",
        nodes=data["nodes"],
        links=data["links"],
        itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#aaa"),
        linestyle_opt=opts.LineStyleOpts(color="source", curve=0.5, opacity=0.5),
        tooltip_opts=opts.TooltipOpts(trigger_on="mousemove"),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Sankey Diagram"))
)

# ---------------------------------------------- Sankey - Sankey_with_level_setting ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Sankey

with open("../assets/echarts/sankey/product.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c3 = (
    Sankey()
    .add(
        "sankey",
        nodes=j["nodes"],
        links=j["links"],
        pos_top="10%",
        levels=[
            opts.SankeyLevelsOpts(
                depth=0,
                itemstyle_opts=opts.ItemStyleOpts(color="#fbb4ae"),
                linestyle_opts=opts.LineStyleOpts(color="source", opacity=0.6),
            ),
            opts.SankeyLevelsOpts(
                depth=1,
                itemstyle_opts=opts.ItemStyleOpts(color="#b3cde3"),
                linestyle_opts=opts.LineStyleOpts(color="source", opacity=0.6),
            ),
            opts.SankeyLevelsOpts(
                depth=2,
                itemstyle_opts=opts.ItemStyleOpts(color="#ccebc5"),
                linestyle_opts=opts.LineStyleOpts(color="source", opacity=0.6),
            ),
            opts.SankeyLevelsOpts(
                depth=3,
                itemstyle_opts=opts.ItemStyleOpts(color="#decbe4"),
                linestyle_opts=opts.LineStyleOpts(color="source", opacity=0.6),
            ),
        ],
        linestyle_opt=opts.LineStyleOpts(curve=0.5),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Sankey-Level Settings"),
        tooltip_opts=opts.TooltipOpts(trigger="item", trigger_on="mousemove"),
    )
)

# ---------------------------------------------- Sankey - Sankey_vertical ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Sankey

colors = [
    "#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061",
]
nodes = [
    {"name": "a"},
    {"name": "b"},
    {"name": "a1"},
    {"name": "b1"},
    {"name": "c"},
    {"name": "e"},
]
links = [
    {"source": "a", "target": "a1", "value": 5},
    {"source": "e", "target": "b", "value": 3},
    {"source": "a", "target": "b1", "value": 3},
    {"source": "b1", "target": "a1", "value": 1},
    {"source": "b1", "target": "c", "value": 2},
    {"source": "b", "target": "c", "value": 1},
]
c4 = (
    Sankey()
    .set_colors(colors)
    .add(
        "sankey",
        nodes=nodes,
        links=links,
        pos_bottom="10%",
        # focus_node_mode="allEdges",
        orient="vertical",
        linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),
        label_opts=opts.LabelOpts(position="top"),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Sankey-Vertical"),
        tooltip_opts=opts.TooltipOpts(trigger="item", trigger_on="mousemove"),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="sankey",
    public=True,
)

swanlab.log(
    {
        "Sankey - Sankey_base": c1,
        "Sankey - Sankey_diagram": c2,
        "Sankey - Sankey_with_level_setting": c3,
        "Sankey - Sankey_vertical": c4
    }
)
