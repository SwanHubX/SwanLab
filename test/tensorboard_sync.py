from tensorboardX import SummaryWriter
# from torch.utils.tensorboard import SummaryWriter
import numpy as np
import swanlab

swanlab.sync_tensorboardX()
# swanlab.sync_tensorboard_torch()

writer = SummaryWriter('runs/example')

for i in range(100):
    scalar_value = np.random.rand()
    writer.add_scalar('random_scalar', scalar_value, i)

writer.close()