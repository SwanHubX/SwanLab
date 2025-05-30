"""
// @author: ComPleHN
// @file: scatter3d.vue
// @time: 2025/5/29 10:35
// @description: 本文件是对于echarts的 3D散点图 图表的测试
"""
# ---------------------------------------------- Scatter3d - Scatter3d ----------------------------------------------
import asyncio
from aiohttp import TCPConnector, ClientSession

import pyecharts.options as opts
from pyecharts.charts import Scatter3D

"""
Gallery 使用 pyecharts 1.1.0
参考地址: https://echarts.apache.org/examples/editor.html?c=scatter3d&gl=1&theme=dark

目前无法实现的功能:

1、暂时无法对 Grid3D 设置 轴线和轴坐标的 style (非白色背景下有问题)
"""


async def get_json_data(url: str) -> dict:
    async with ClientSession(connector=TCPConnector(ssl=False)) as session:
        async with session.get(url=url) as response:
            return await response.json()


# 获取官方的数据
data = asyncio.run(
    get_json_data(
        url="https://echarts.apache.org/examples/data/asset/data/nutrients.json"
    )
)

# 列名映射
field_indices = {
    "calcium": 3,
    "calories": 12,
    "carbohydrate": 8,
    "fat": 10,
    "fiber": 5,
    "group": 1,
    "id": 16,
    "monounsat": 14,
    "name": 0,
    "polyunsat": 15,
    "potassium": 7,
    "protein": 2,
    "saturated": 13,
    "sodium": 4,
    "sugars": 9,
    "vitaminc": 6,
    "water": 11,
}

# 配置 config
config_xAxis3D = "protein"
config_yAxis3D = "fiber"
config_zAxis3D = "sodium"
config_color = "fiber"
config_symbolSize = "vitaminc"

# 构造数据
data = [
    [
        item[field_indices[config_xAxis3D]],
        item[field_indices[config_yAxis3D]],
        item[field_indices[config_zAxis3D]],
        item[field_indices[config_color]],
        item[field_indices[config_symbolSize]],
        index,
    ]
    for index, item in enumerate(data)
]

c1 = (
    Scatter3D()  # bg_color="black"
    .add(
        series_name="",
        data=data,
        xaxis3d_opts=opts.Axis3DOpts(
            name=config_xAxis3D,
            type_="value",
            # textstyle_opts=opts.TextStyleOpts(color="#fff"),
        ),
        yaxis3d_opts=opts.Axis3DOpts(
            name=config_yAxis3D,
            type_="value",
            # textstyle_opts=opts.TextStyleOpts(color="#fff"),
        ),
        zaxis3d_opts=opts.Axis3DOpts(
            name=config_zAxis3D,
            type_="value",
            # textstyle_opts=opts.TextStyleOpts(color="#fff"),
        ),
        grid3d_opts=opts.Grid3DOpts(width=100, height=100, depth=100),
    )
    .set_global_opts(
        visualmap_opts=[
            opts.VisualMapOpts(
                type_="color",
                is_calculable=True,
                dimension=3,
                pos_top="10",
                max_=79 / 2,
                range_color=[
                    "#1710c0",
                    "#0b9df0",
                    "#00fea8",
                    "#00ff0d",
                    "#f5f811",
                    "#f09a09",
                    "#fe0300",
                ],
            ),
            opts.VisualMapOpts(
                type_="size",
                is_calculable=True,
                dimension=4,
                pos_bottom="10",
                max_=2.4 / 2,
                range_size=[10, 40],
            ),
        ]
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="scatter3d",
    public=True,
)

swanlab.log(
    {
        "Scatter3d - Scatter3d": c1
    }
)
