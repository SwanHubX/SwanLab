import swanlab

def sync_tensorboardX():
    """
    同步tensorboardX到swanlab
    
    from tensorboardX import SummaryWriter
    import numpy as np

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
    """
    try:
        from tensorboardX import SummaryWriter
    except ImportError:
        raise ImportError("please install tensorboardX first, command: `pip install tensorboardX`")
    
    original_init = SummaryWriter.__init__
    original_add_scalar = SummaryWriter.add_scalar
    original_close = SummaryWriter.close
    
    def patched_init(self, *args, **kwargs):
        logdir = args[0]
            
        tb_config = {
            'tensorboard_logdir': logdir,
        }
        
        if swanlab.data.get_run() is None:
            swanlab.init(config=tb_config)
        else:
            swanlab.config.update(tb_config)
        
        return original_init(self, *args, **kwargs)
    
    def patched_add_scalar(self, tag, scalar_value, global_step=None):
        data = {tag: scalar_value}
        swanlab.log(data=data, step=global_step)
        
        return original_add_scalar(self, tag, scalar_value, global_step)
            
    def patched_close(self):
        # 调用原始的close方法
        original_close(self)
        # 关闭SwanLab记录器
        swanlab.finish()
    
    # 应用monkey patch
    SummaryWriter.__init__ = patched_init
    SummaryWriter.add_scalar = patched_add_scalar
    SummaryWriter.close = patched_close


def sync_tensorboard_torch():
    """
    同步torch自带的tensorboard到swanlab
    
    from torch.utils.tensorboard import SummaryWriter
    import numpy as np

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
    """
    try:
        from torch.utils.tensorboard import SummaryWriter
    except ImportError:
        raise ImportError("please install torch first, command: `pip install torch`")

    original_init = SummaryWriter.__init__
    original_add_scalar = SummaryWriter.add_scalar
    original_close = SummaryWriter.close
    
    def patched_init(self, *args, **kwargs):
        logdir = args[0]
            
        tb_config = {
            'tensorboard_logdir': logdir,
        }
        
        if swanlab.data.get_run() is None:
            swanlab.init(config=tb_config)
        else:
            swanlab.config.update(tb_config)
        
        return original_init(self, *args, **kwargs)
    
    def patched_add_scalar(self, tag, scalar_value, global_step=None):
        data = {tag: scalar_value}
        swanlab.log(data=data, step=global_step)
        
        return original_add_scalar(self, tag, scalar_value, global_step)
            
    def patched_close(self):
        # 调用原始的close方法
        original_close(self)
        # 关闭SwanLab记录器
        swanlab.finish()
    
    # 应用monkey patch
    SummaryWriter.__init__ = patched_init
    SummaryWriter.add_scalar = patched_add_scalar
    SummaryWriter.close = patched_close

    