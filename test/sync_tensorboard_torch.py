"""
@author: cunyue
@file: sync_tensorboard_torch.py
@time: 2025/1/20 17:15
@description: 测试同步tensorboard torch
"""

import numpy as np
from torch.utils.tensorboard import SummaryWriter

import swanlab

# 只同步标量数据（排除文本、图像等）
swanlab.sync_tensorboard_torch(types=['scalar', 'scalars'])
writer = SummaryWriter('runs/example')

# 这些不会被同步（因为 types 过滤）
writer.add_image('random_image', np.random.randint(0, 255, (3, 100, 100)), global_step=20)
writer.add_text('random_text', 'hello', global_step=10)

for i in range(100):
    scalar_value = np.random.rand()
    # 这些会被同步（标量数据）
    writer.add_scalar('random_scalar', scalar_value, i)
    writer.add_scalars('random_scalars', {'scalar1': scalar_value, 'scalar2': scalar_value * 2}, i)

writer.close()