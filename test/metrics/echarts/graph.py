"""
// @author: ComPleHN
// @file: graph.vue
// @time: 2025/5/27 16:04
// @description: 本文件是对于echarts的 关系图 图表测试,本文件有注意事项请看本目录下 README.md
"""
# ---------------------------------------------- Graph - Graph_base ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Graph

nodes = [
    {"name": "结点1", "symbolSize": 10},
    {"name": "结点2", "symbolSize": 20},
    {"name": "结点3", "symbolSize": 30},
    {"name": "结点4", "symbolSize": 40},
    {"name": "结点5", "symbolSize": 50},
    {"name": "结点6", "symbolSize": 40},
    {"name": "结点7", "symbolSize": 30},
    {"name": "结点8", "symbolSize": 20},
]
links = []
for i in nodes:
    for j in nodes:
        links.append({"source": i.get("name"), "target": j.get("name")})
c1 = (
    Graph()
    .add("", nodes, links, repulsion=8000)
    .set_global_opts(title_opts=opts.TitleOpts(title="Graph-基本示例"))
)

# ---------------------------------------------- Graph - Graph_with_edge_options ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Graph

nodes_data = [
    opts.GraphNode(name="结点1", symbol_size=10),
    opts.GraphNode(name="结点2", symbol_size=20),
    opts.GraphNode(name="结点3", symbol_size=30),
    opts.GraphNode(name="结点4", symbol_size=40),
    opts.GraphNode(name="结点5", symbol_size=50),
    opts.GraphNode(name="结点6", symbol_size=60),
]
links_data = [
    opts.GraphLink(source="结点1", target="结点2", value=2),
    opts.GraphLink(source="结点2", target="结点3", value=3),
    opts.GraphLink(source="结点3", target="结点4", value=4),
    opts.GraphLink(source="结点4", target="结点5", value=5),
    opts.GraphLink(source="结点5", target="结点6", value=6),
    opts.GraphLink(source="结点6", target="结点1", value=7),
]
c2 = (
    Graph()
    .add(
        "",
        nodes_data,
        links_data,
        repulsion=4000,
        edge_label=opts.LabelOpts(
            is_show=True, position="middle", formatter="{b} 的数据 {c}"
        ),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Graph-GraphNode-GraphLink-WithEdgeLabel")
    )
)

# ---------------------------------------------- Graph - Graph_with_options ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Graph

nodes = [
    opts.GraphNode(name="结点1", symbol_size=10),
    opts.GraphNode(name="结点2", symbol_size=20),
    opts.GraphNode(name="结点3", symbol_size=30),
    opts.GraphNode(name="结点4", symbol_size=40),
    opts.GraphNode(name="结点5", symbol_size=50),
]
links = [
    opts.GraphLink(source="结点1", target="结点2"),
    opts.GraphLink(source="结点2", target="结点3"),
    opts.GraphLink(source="结点3", target="结点4"),
    opts.GraphLink(source="结点4", target="结点5"),
    opts.GraphLink(source="结点5", target="结点1"),
]
c3 = (
    Graph()
    .add("", nodes, links, repulsion=4000)
    .set_global_opts(title_opts=opts.TitleOpts(title="Graph-GraphNode-GraphLink"))
)

# ---------------------------------------------- Graph - Graph_weibo ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Graph

with open("../assets/echarts/graph/weibo.json", "r", encoding="utf-8") as f:
    j = json.load(f)
    nodes, links, categories, cont, mid, userl = j
c4 = (
    Graph()
    .add(
        "",
        nodes,
        links,
        categories,
        repulsion=50,
        linestyle_opts=opts.LineStyleOpts(curve=0.2),
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        legend_opts=opts.LegendOpts(is_show=False),
        title_opts=opts.TitleOpts(title="Graph-微博转发关系图"),
    )
)

# ---------------------------------------------- Graph - Npm_dependencies ---------------------------------------------- 
import asyncio
from aiohttp import TCPConnector, ClientSession

import pyecharts.options as opts
from pyecharts.charts import Graph

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=graph-npm

目前无法实现的功能:

1、暂无
"""


async def get_json_data(url: str) -> dict:
    async with ClientSession(connector=TCPConnector(ssl=False)) as session:
        async with session.get(url=url) as response:
            return await response.json()


# 获取官方的数据
data = asyncio.run(
    get_json_data(
        url="https://echarts.apache.org/examples/data/asset/data/npmdepgraph.min10.json"
    )
)

nodes = [
    {
        "x": node["x"],
        "y": node["y"],
        "id": node["id"],
        "name": node["label"],
        "symbolSize": node["size"],
        "itemStyle": {"normal": {"color": node["color"]}},
    }
    for node in data["nodes"]
]

edges = [
    {"source": edge["sourceID"], "target": edge["targetID"]} for edge in data["edges"]
]


c5 = (
    Graph()
    .add(
        series_name="",
        nodes=nodes,
        links=edges,
        layout="none",
        is_roam=True,
        is_focusnode=True,
        label_opts=opts.LabelOpts(is_show=False),
        linestyle_opts=opts.LineStyleOpts(width=0.5, curve=0.3, opacity=0.7),
    )
    .set_global_opts(title_opts=opts.TitleOpts(title="NPM Dependencies"))
)

# ---------------------------------------------- Graph - Graph_les_miserables ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Graph


with open("../assets/echarts/graph/les-miserables.json", "r", encoding="utf-8") as f:
    j = json.load(f)
    nodes = j["nodes"]
    links = j["links"]
    categories = j["categories"]

c6 = (
    Graph(init_opts=opts.InitOpts(width="1000px", height="600px"))
    .add(
        "",
        nodes=nodes,
        links=links,
        categories=categories,
        layout="circular",
        is_rotate_label=True,
        linestyle_opts=opts.LineStyleOpts(color="source", curve=0.3),
        label_opts=opts.LabelOpts(position="right"),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Graph-Les Miserables"),
        legend_opts=opts.LegendOpts(orient="vertical", pos_left="2%", pos_top="20%"),
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="graph",
    public=True,
)

swanlab.log(
    {
        "Graph - Graph_base": c1,
        "Graph - Graph_with_edge_options": c2,
        "Graph - Graph_with_options": c3,
        "Graph - Graph_weibo": c4,
        "Graph - Npm_dependencies": c5,
        "Graph - Graph_les_miserables": c6
    }
)

