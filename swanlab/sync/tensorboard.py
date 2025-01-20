import swanlab


def sync_tensorboardX():
    """
    同步tensorboardX到swanlab

    from tensorboardX import SummaryWriter
    import numpy as np
    import swanlab

    swanlab.sync_tensorboardX()
    writer = SummaryWriter('runs/example')

    for i in range(100):
        scalar_value = np.random.rand()
        writer.add_scalar('random_scalar', scalar_value, i)

    writer.close()
    """
    try:
        from tensorboardX import SummaryWriter
    except ImportError:
        raise ImportError("please install tensorboardX first, command: `pip install tensorboardX`")

    original_init = SummaryWriter.__init__
    original_add_scalar = SummaryWriter.add_scalar
    original_add_image = SummaryWriter.add_image
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

    def patched_add_image(self, tag, img_tensor, global_step=None, walltime=None, dataformats='CHW'):
        import numpy as np
        
        img_tensor_swanlab = img_tensor.copy()
        # Convert to numpy array if it's a tensor
        if hasattr(img_tensor, 'cpu'):
            img_tensor_swanlab = img_tensor.cpu()
        if hasattr(img_tensor, 'numpy'):
            img_tensor_swanlab = img_tensor.numpy()
            
        # Handle different input formats
        if dataformats == 'CHW':
            # Convert CHW to HWC for swanlab
            img_tensor_swanlab = np.transpose(img_tensor_swanlab, (1, 2, 0))
        elif dataformats == 'NCHW':
            # Take first image if batch dimension exists and convert to HWC
            img_tensor_swanlab = np.transpose(img_tensor_swanlab[0], (1, 2, 0))
        elif dataformats == 'HW':
            # Add channel dimension for grayscale
            img_tensor_swanlab = np.expand_dims(img_tensor_swanlab, axis=-1)
        elif dataformats == 'HWC':
            # Already in correct format
            pass
            
        data = {tag: swanlab.Image(img_tensor_swanlab)}
        swanlab.log(data=data, step=global_step)

        return original_add_image(self, tag, img_tensor, global_step, walltime, dataformats)

    def patched_close(self):
        # 调用原始的close方法
        original_close(self)
        # 关闭SwanLab记录器
        swanlab.finish()

    # 应用monkey patch
    SummaryWriter.__init__ = patched_init
    SummaryWriter.add_scalar = patched_add_scalar
    SummaryWriter.add_image = patched_add_image
    SummaryWriter.close = patched_close


def sync_tensorboard_torch():
    """
    同步torch自带的tensorboard到swanlab

    from torch.utils.tensorboard import SummaryWriter
    import numpy as np
    import swanlab

    swanlab.sync_tensorboard_torch()
    writer = SummaryWriter('runs/example')

    for i in range(100):
        scalar_value = np.random.rand()
        writer.add_scalar('random_scalar', scalar_value, i)

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
