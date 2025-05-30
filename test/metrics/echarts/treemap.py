"""
// @author: ComPleHN
// @file: treemap.vue
// @time: 2025/5/29 12:35
// @description: 本文件是对于echarts的 矩形树图 图表的测试,本文件有注意事项请看本目录下 README.md
"""
# ---------------------------------------------- Treemap - Echarts_option_query ----------------------------------------------
import re
import asyncio
from aiohttp import TCPConnector, ClientSession

import pyecharts.options as opts
from pyecharts.charts import TreeMap

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=treemap-drill-down

目前无法实现的功能:

1、层级的样式配置
"""


async def get_json_data(url: str) -> dict:
    async with ClientSession(connector=TCPConnector(ssl=False)) as session:
        async with session.get(url=url) as response:
            return await response.json()


# 获取官方的数据
data = asyncio.run(
    get_json_data(
        url="https://echarts.apache.org/examples/data/asset/data/"
        "ec-option-doc-statistics-201604.json"
    )
)

tree_map_data: dict = {"children": []}


def convert(source, target, base_path: str):
    for key in source:
        if base_path != "":
            path = base_path + "." + key
        else:
            path = key
        if re.match(r"/^\$/", key):
            pass
        else:
            child = {"name": path, "children": []}
            target["children"].append(child)
            if isinstance(source[key], dict):
                convert(source[key], child, path)
            else:
                target["value"] = source["$count"]


convert(source=data, target=tree_map_data, base_path="")


c1 = (
    TreeMap(init_opts=opts.InitOpts(width="1200px", height="720px"))
    .add(
        series_name="option",
        data=tree_map_data["children"],
        visual_min=300,
        leaf_depth=1,
        # 标签居中为 position = "inside"
        label_opts=opts.LabelOpts(position="inside"),
    )
    .set_global_opts(
        legend_opts=opts.LegendOpts(is_show=False),
        title_opts=opts.TitleOpts(
            title="Echarts 配置项查询分布", subtitle="2016/04", pos_left="leafDepth"
        ),
    )
)

# ---------------------------------------------- Treemap - Treemap_levels ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import TreeMap


with open("../assets/echarts/treemap/treemap.json", "r", encoding="utf-8") as f:
    data = json.load(f)
c2 = (
    TreeMap()
    .add(
        series_name="演示数据",
        data=data,
        levels=[
            opts.TreeMapLevelsOpts(
                treemap_itemstyle_opts=opts.TreeMapItemStyleOpts(
                    border_color="#555", border_width=4, gap_width=4
                )
            ),
            opts.TreeMapLevelsOpts(
                color_saturation=[0.3, 0.6],
                treemap_itemstyle_opts=opts.TreeMapItemStyleOpts(
                    border_color_saturation=0.7, gap_width=2, border_width=2
                ),
            ),
            opts.TreeMapLevelsOpts(
                color_saturation=[0.3, 0.5],
                treemap_itemstyle_opts=opts.TreeMapItemStyleOpts(
                    border_color_saturation=0.6, gap_width=1
                ),
            ),
            opts.TreeMapLevelsOpts(color_saturation=[0.3, 0.5]),
        ],
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="TreeMap-Levels-配置"))
)

# ---------------------------------------------- Treemap - Treemap_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import TreeMap

data = [
    {"value": 40, "name": "我是A"},
    {
        "value": 180,
        "name": "我是B",
        "children": [
            {
                "value": 76,
                "name": "我是B.children",
                "children": [
                    {"value": 12, "name": "我是B.children.a"},
                    {"value": 28, "name": "我是B.children.b"},
                    {"value": 20, "name": "我是B.children.c"},
                    {"value": 16, "name": "我是B.children.d"},
                ],
            }
        ],
    },
]

c3 = (
    TreeMap()
    .add("演示数据", data)
    .set_global_opts(title_opts=opts.TitleOpts(title="TreeMap-基本示例"))
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="treemap",
    public=True,
)

swanlab.log(
    {
        "Treemap - Echarts_option_query": c1,
        "Treemap - Treemap_levels": c2,
        "Treemap - Treemap_base": c3
    }
)
