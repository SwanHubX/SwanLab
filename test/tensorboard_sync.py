# from tensorboardX import SummaryWriter
from torch.utils.tensorboard import SummaryWriter
import swanlab
import numpy as np

swanlab.sync_tensorboard_torch()

# 创建一个SummaryWriter对象，指定日志目录
writer = SummaryWriter('runs/example')

# 生成一些示例数据
for i in range(100):
    # 生成一个随机的标量值
    scalar_value = np.random.rand()
    
    # 记录标量值
    writer.add_scalar('random_scalar', scalar_value, i)

# 关闭SummaryWriter
writer.close()