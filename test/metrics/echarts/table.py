"""
// @author: ComPleHN
// @file: table.vue
// @time: 2025/5/29 10:46
// @description: 本文件是对于echarts的 表格组件 图表的测试
"""
# ---------------------------------------------- Table - Table_base ----------------------------------------------
from pyecharts.components import Table
from pyecharts.options import ComponentTitleOpts


table = Table()

headers = ["City name", "Area", "Population", "Annual Rainfall"]
rows = [
    ["Brisbane", 5905, 1857594, 1146.4],
    ["Adelaide", 1295, 1158259, 600.5],
    ["Darwin", 112, 120900, 1714.7],
    ["Hobart", 1357, 205556, 619.5],
    ["Sydney", 2058, 4336374, 1214.8],
    ["Melbourne", 1566, 3806092, 646.9],
    ["Perth", 5386, 1554769, 869.4],
]
table.add(headers, rows)
table.set_global_opts(
    title_opts=ComponentTitleOpts(title="Table-基本示例", subtitle="我是副标题支持换行哦")
)

c1 = table

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="table",
    public=True,
)

swanlab.log(
    {
        "Table - Table_base": c1
    }
)
