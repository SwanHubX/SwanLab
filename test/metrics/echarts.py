"""
@author: cunyue
@file: echarts.py
@time: 2025/5/21 14:33
@description: echarts 模块的测试
"""

import swanlab

swanlab.init(project="echarts", public=True)


bar = swanlab.echarts.Bar().add_xaxis(["A", "B", "C"]).add_yaxis("数据", [1, 2, 3])

swanlab.log({"bar": bar})
