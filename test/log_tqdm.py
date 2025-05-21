"""
@author: cunyue
@file: proxy_tqdm.py
@time: 2025/5/17 21:23
@description: 测试tqdm的进度条在 swanlab 终端代理下的表现
"""

import time

from tqdm import tqdm

import swanlab

swanlab.init(project='log-tqdm')

# 在循环中显示进度条
for i in tqdm(range(100)):
    time.sleep(0.01)
    if i % 10 == 0:
        swanlab.log({"progress": i}, print_to_console=True)
