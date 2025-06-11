"""
// @author: ComPleHN
// @file: table.vue
// @time: 2025/5/29 10:46
// @description: 本文件是对于echarts的 表格组件 图表的测试
"""
# ---------------------------------------------- Table - Table_base ----------------------------------------------
from swanlab.data.modules.custom_charts.table import Table
from pyecharts.options import ComponentTitleOpts


table = Table()

headers = ["NO", "机构", "数量"]
rows = [
    [1, "个人", 10],
    [2, "西安电子科技大学", 259 ],
    [3, "西安邮电大学", 123 ],
    [4, "北京大学", 300 ],
    [5, "清华大学", 290],
    [6, "helloworld", 1145]
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
