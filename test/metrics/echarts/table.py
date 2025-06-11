"""
// @author: ComPleHN
// @file: table.vue
// @time: 2025/5/29 10:46
// @description: 本文件是对于echarts的 表格组件 图表的测试
"""

# ---------------------------------------------- Table - Table_base ----------------------------------------------
from swanlab import echarts


table = echarts.Table()

headers = ["NO", "Product", "Count"]
rows = [
    [2, "A", 259],
    [3, "B", 123],
    [4, "C", 300],
    [5, "D", 290],
    [6, "E", 1145],
]
table.add(headers, rows)

c1 = table

import swanlab

swanlab.init(
    project="echarts-test",
    experiment_name="table",
    public=True,
)

swanlab.log({"Table - Table_base": c1})
