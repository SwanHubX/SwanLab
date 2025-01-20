"""
@author: cunyue
@file: sync_tensorboardX.py
@time: 2025/1/20 17:15
@description: 测试同步tensorboardX
"""

import numpy as np
from tensorboardX import SummaryWriter

import swanlab

swanlab.sync_tensorboardX()

writer = SummaryWriter('runs/example')

for i in range(100):
    scalar_value = np.random.rand()
    writer.add_scalar('random_scalar', scalar_value, i)

writer.close()
