"""
// @author: ComPleHN
// @file: tree.vue
// @time: 2025/5/29 12:20
// @description: 本文件是对于echarts的 树状图 图表的测试,本文件有注意事项请看本目录下 README.md
"""
# ---------------------------------------------- Tree - Tree_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Tree


data = [
    {
        "children": [
            {"name": "B"},
            {
                "children": [{"children": [{"name": "I"}], "name": "E"}, {"name": "F"}],
                "name": "C",
            },
            {
                "children": [
                    {"children": [{"name": "J"}, {"name": "K"}], "name": "G"},
                    {"name": "H"},
                ],
                "name": "D",
            },
        ],
        "name": "A",
    }
]
c1 = (
    Tree()
    .add("", data)
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-基本示例"))
)

# ---------------------------------------------- Tree - Tree_bottom_top ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Tree

with open("../assets/echarts/tree/flare.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c2 = (
    Tree()
    .add(
        "",
        [j],
        collapse_interval=2,
        orient="BT",
        label_opts=opts.LabelOpts(
            position="top",
            horizontal_align="right",
            vertical_align="middle",
            rotate=-90,
        ),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-下上方向"))
)

# ---------------------------------------------- Tree - Tree_right_left ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Tree

with open("../assets/echarts/tree/flare.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c3 = (
    Tree()
    .add("", [j], collapse_interval=2, orient="RL")
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-右左方向"))
)

# ---------------------------------------------- Tree - Tree_left_right ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Tree


with open("../assets/echarts/tree/flare.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c4 = (
    Tree()
    .add("", [j], collapse_interval=2)
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-左右方向"))
)

# ---------------------------------------------- Tree - Tree_top_bottom ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Tree


with open("../assets/echarts/tree/flare.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c5 = (
    Tree()
    .add(
        "",
        [j],
        collapse_interval=2,
        orient="TB",
        label_opts=opts.LabelOpts(
            position="top",
            horizontal_align="right",
            vertical_align="middle",
            rotate=-90,
        ),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-上下方向"))
)

# ---------------------------------------------- Tree - Radial_tree ----------------------------------------------
import asyncio
from aiohttp import TCPConnector, ClientSession

import pyecharts.options as opts
from pyecharts.charts import Tree

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=tree-radial

目前无法实现的功能:

1、
"""


async def get_json_data(url: str) -> dict:
    async with ClientSession(connector=TCPConnector(ssl=False)) as session:
        async with session.get(url=url) as response:
            return await response.json()


# 获取官方的数据
data = asyncio.run(
    get_json_data(url="https://echarts.apache.org/examples/data/asset/data/flare.json")
)

c6 = (
    Tree()
    .add(
        series_name="",
        data=[data],
        pos_top="18%",
        pos_bottom="14%",
        layout="radial",
        symbol="emptyCircle",
        symbol_size=7,
    )
    .set_global_opts(
        tooltip_opts=opts.TooltipOpts(trigger="item", trigger_on="mousemove")
    )
)

# ---------------------------------------------- Tree - Tree_layout ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Tree

with open("../assets/echarts/tree/flare.json", "r", encoding="utf-8") as f:
    j = json.load(f)
c7 = (
    Tree()
    .add("", [j], collapse_interval=2, layout="radial")
    .set_global_opts(title_opts=opts.TitleOpts(title="Tree-Layout"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="tree",
    public=True,
)

swanlab.log(
    {
        "Tree - Tree_base": c1,
        "Tree - Tree_bottom_top": c2,
        "Tree - Tree_right_left": c3,
        "Tree - Tree_left_right": c4,
        "Tree - Tree_top_bottom": c5,
        "Tree - Radial_tree": c6,
        "Tree - Tree_layout": c7
    }
)
