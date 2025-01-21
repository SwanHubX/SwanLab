"""
@author: cunyue
@file: sync_tensorboard_torch.py
@time: 2025/1/20 17:15
@description: 测试同步tensorboard torch
"""

import numpy as np
from torch.utils.tensorboard import SummaryWriter

import swanlab

swanlab.sync_tensorboard_torch()
writer = SummaryWriter('runs/example')

writer.add_image('random_image', np.random.randint(0, 255, (3, 100, 100)), global_step=20)

for i in range(100):
    scalar_value = np.random.rand()
    writer.add_scalar('random_scalar', scalar_value, i)

writer.close()
