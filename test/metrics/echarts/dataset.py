"""
// @author: ComPleHN
// @file: dataset.vue
// @time: 2025/5/27 15:03
// @description: 本文件是对于echarts的 数据集 图表的测试,本文件有注意事项请看本目录下 README.md
"""

# ---------------------------------------------- Dataset - Dataset_pie ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Pie

c1 = (
    Pie()
    .add_dataset(
        source=[
            ["product", "2012", "2013", "2014", "2015", "2016", "2017"],
            ["Matcha Latte", 41.1, 30.4, 65.1, 53.3, 83.8, 98.7],
            ["Milk Tea", 86.5, 92.1, 85.7, 83.1, 73.4, 55.1],
            ["Cheese Cocoa", 24.1, 67.2, 79.5, 86.4, 65.2, 82.5],
            ["Walnut Brownie", 55.2, 67.1, 69.2, 72.4, 53.9, 39.1],
        ]
    )
    .add(
        series_name="Matcha Latte",
        data_pair=[],
        radius=60,
        center=["25%", "30%"],
        encode={"itemName": "product", "value": "2012"},
    )
    .add(
        series_name="Milk Tea",
        data_pair=[],
        radius=60,
        center=["75%", "30%"],
        encode={"itemName": "product", "value": "2013"},
    )
    .add(
        series_name="Cheese Cocoa",
        data_pair=[],
        radius=60,
        center=["25%", "75%"],
        encode={"itemName": "product", "value": "2014"},
    )
    .add(
        series_name="Walnut Brownie",
        data_pair=[],
        radius=60,
        center=["75%", "75%"],
        encode={"itemName": "product", "value": "2015"},
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Dataset simple pie example"),
        legend_opts=opts.LegendOpts(pos_left="30%", pos_top="2%"),
    )
)

# ---------------------------------------------- Dataset - Dataset_bar_0 ----------------------------------------------
from pyecharts import options as opts
from pyecharts.charts import Bar

c2 = (
    Bar()
    .add_dataset(
        source=[
            ["product", "2015", "2016", "2017"],
            ["Matcha Latte", 43.3, 85.8, 93.7],
            ["Milk Tea", 83.1, 73.4, 55.1],
            ["Cheese Cocoa", 86.4, 65.2, 82.5],
            ["Walnut Brownie", 72.4, 53.9, 39.1],
        ]
    )
    .add_yaxis(series_name="2015", y_axis=[])
    .add_yaxis(series_name="2016", y_axis=[])
    .add_yaxis(series_name="2017", y_axis=[])
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Dataset simple bar example"),
        xaxis_opts=opts.AxisOpts(type_="category"),
    )
)

# ---------------------------------------------- Dataset - Dataset_bar_1 ---------------------------------------------- 
from pyecharts import options as opts
from pyecharts.charts import Bar

c3 = (
    Bar()
    .add_dataset(
        source=[
            ["score", "amount", "product"],
            [89.3, 58212, "Matcha Latte"],
            [57.1, 78254, "Milk Tea"],
            [74.4, 41032, "Cheese Cocoa"],
            [50.1, 12755, "Cheese Brownie"],
            [89.7, 20145, "Matcha Cocoa"],
            [68.1, 79146, "Tea"],
            [19.6, 91852, "Orange Juice"],
            [10.6, 101852, "Lemon Juice"],
            [32.7, 20112, "Walnut Brownie"],
        ]
    )
    .add_yaxis(
        series_name="",
        y_axis=[],
        encode={"x": "amount", "y": "product"},
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        title_opts=opts.TitleOpts(title="Dataset normal bar example"),
        xaxis_opts=opts.AxisOpts(name="amount"),
        yaxis_opts=opts.AxisOpts(type_="category"),
        visualmap_opts=opts.VisualMapOpts(
            orient="horizontal",
            pos_left="center",
            min_=10,
            max_=100,
            range_text=["High Score", "Low Score"],
            dimension=0,
            range_color=["#D7DA8B", "#E15457"],
        ),
    )
)

# ---------------------------------------------- Dataset - Dataset_professional_scatter ----------------------------------------------
import json

from pyecharts import options as opts
from pyecharts.charts import Grid, Scatter

with open("../assets/echarts/dataset/life-expectancy-table.json", "r", encoding="utf-8") as f:
    j = json.load(f)

l1_1 = (
    Scatter()
    .add_dataset(
        dimensions=[
            "Income",
            "Life Expectancy",
            "Population",
            "Country",
            {"name": "Year", "type": "ordinal"},
        ],
        source=j,
    )
    .add_yaxis(
        series_name="",
        y_axis=[],
        symbol_size=2.5,
        xaxis_index=0,
        yaxis_index=0,
        encode={"x": "Income", "y": "Life Expectancy", "tooltip": [0, 1, 2, 3, 4]},
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value",
            grid_index=0,
            name="Income",
            axislabel_opts=opts.LabelOpts(rotate=50, interval=0),
        ),
        yaxis_opts=opts.AxisOpts(type_="value", grid_index=0, name="Life Expectancy"),
        title_opts=opts.TitleOpts(title="Encode and Matrix"),
    )
)

l1_2 = (
    Scatter()
    .add_dataset()
    .add_yaxis(
        series_name="",
        y_axis=[],
        symbol_size=2.5,
        xaxis_index=1,
        yaxis_index=1,
        encode={"x": "Country", "y": "Income", "tooltip": [0, 1, 2, 3, 4]},
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="category",
            grid_index=1,
            name="Country",
            boundary_gap=False,
            axislabel_opts=opts.LabelOpts(rotate=50, interval=0),
        ),
        yaxis_opts=opts.AxisOpts(type_="value", grid_index=1, name="Income"),
    )
)

l2_1 = (
    Scatter()
    .add_dataset()
    .add_yaxis(
        series_name="",
        y_axis=[],
        symbol_size=2.5,
        xaxis_index=2,
        yaxis_index=2,
        encode={"x": "Income", "y": "Population", "tooltip": [0, 1, 2, 3, 4]},
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value",
            grid_index=2,
            name="Income",
            axislabel_opts=opts.LabelOpts(rotate=50, interval=0),
        ),
        yaxis_opts=opts.AxisOpts(type_="value", grid_index=2, name="Population"),
    )
)

l2_2 = (
    Scatter()
    .add_dataset()
    .add_yaxis(
        series_name="",
        y_axis=[],
        symbol_size=2.5,
        xaxis_index=3,
        yaxis_index=3,
        encode={"x": "Life Expectancy", "y": "Population", "tooltip": [0, 1, 2, 3, 4]},
        label_opts=opts.LabelOpts(is_show=False),
    )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value",
            grid_index=3,
            name="Life Expectancy",
            axislabel_opts=opts.LabelOpts(rotate=50, interval=0),
        ),
        yaxis_opts=opts.AxisOpts(type_="value", grid_index=3, name="Population"),
    )
)

c4 = (
    Grid(init_opts=opts.InitOpts(width="1200px", height="960px"))
    .add(
        chart=l1_1,
        grid_opts=opts.GridOpts(pos_right="57%", pos_bottom="57%"),
        grid_index=0,
    )
    .add(
        chart=l1_2,
        grid_opts=opts.GridOpts(pos_left="57%", pos_bottom="57%"),
        grid_index=1,
    )
    .add(
        chart=l2_1,
        grid_opts=opts.GridOpts(pos_right="57%", pos_top="57%"),
        grid_index=2,
    )
    .add(
        chart=l2_2, grid_opts=opts.GridOpts(pos_left="57%", pos_top="57%"), grid_index=3
    )
)

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="dataset",
    public=True,
)

swanlab.log(
    {
        "Dataset - Dataset_pie": c1,
        "Dataset - Dataset_bar_0": c2,
        "Dataset - Dataset_bar_1": c3,
        "Dataset - Dataset_professional_scatter": c4
    }
)

