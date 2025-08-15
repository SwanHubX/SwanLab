import swanlab

def _extract_args(args, kwargs, param_names):
    """
    从args和kwargs中提取参数值的通用函数
    
    Args:
        args: 位置参数元组
        kwargs: 关键字参数字典
        param_names: 参数名称列表
    
    Returns:
        tuple: 按param_names顺序返回提取的参数值
    """
    values = []
    for i, name in enumerate(param_names):
        if len(args) > i:
            values.append(args[i])
        else:
            values.append(kwargs.get(name, None))
    return tuple(values)


def _create_patched_methods(SummaryWriter, logdir_extractor):
    """
    创建patched方法的工厂函数
    
    Args:
        SummaryWriter: SummaryWriter类
        logdir_extractor: 提取logdir的函数
    
    Returns:
        tuple: (patched_init, patched_add_scalar, patched_add_image, patched_close)
    """
    original_init = SummaryWriter.__init__
    original_add_scalar = SummaryWriter.add_scalar
    original_add_scalars = SummaryWriter.add_scalars
    original_add_image = SummaryWriter.add_image
    original_add_text = SummaryWriter.add_text
    original_close = SummaryWriter.close

    def patched_init(self, *args, **kwargs):
        tb_logdir = logdir_extractor(args, kwargs)
        
        tb_config = {
            'tensorboard_logdir': tb_logdir,
        }

        if swanlab.data.get_run() is None:
            swanlab.init(config=tb_config)
        else:
            swanlab.config.update(tb_config)

        return original_init(self, *args, **kwargs)

    def patched_add_scalar(self, *args, **kwargs):
        tag, scalar_value, global_step = _extract_args(
            args, kwargs, ['tag', 'scalar_value', 'global_step']
        )
        
        data = {tag: scalar_value}
        swanlab.log(data=data, step=int(global_step))
        
        return original_add_scalar(self, *args, **kwargs)

    def patched_add_scalars(self, *args, **kwargs):
        # writer.add_scalars('Loss', {'train': loss_train, 'val': loss_val}, global_step=step)
        tag, scalar_value_dict, global_step = _extract_args(
            args, kwargs, ['tag', 'scalar_value_dict', 'global_step']
        )
        for dict_tag, value in scalar_value_dict.items():
            data = {f"{tag}/{dict_tag}": value}
            swanlab.log(data=data, step=int(global_step))
        return original_add_scalars(self, *args, **kwargs)

    def patched_add_image(self, *args, **kwargs):
        import numpy as np
        
        tag, img_tensor, global_step, dataformats = _extract_args(
            args, kwargs, ['tag', 'img_tensor', 'global_step', 'dataformats']
        )
        dataformats = dataformats or 'CHW'  # 设置默认值
        
        # Convert to numpy array if it's a tensor
        if hasattr(img_tensor, 'cpu'):
            img_tensor = img_tensor.cpu()
        if hasattr(img_tensor, 'numpy'):
            img_tensor = img_tensor.numpy()
            
        # Handle different input formats
        if dataformats == 'CHW':
            # Convert CHW to HWC for swanlab
            img_tensor = np.transpose(img_tensor, (1, 2, 0))
        elif dataformats == 'NCHW':
            # Take first image if batch dimension exists and convert to HWC
            img_tensor = np.transpose(img_tensor, (1, 2, 0))
        elif dataformats == 'HW':
            # Add channel dimension for grayscale
            img_tensor = np.expand_dims(img_tensor, axis=-1)
        elif dataformats == 'HWC':
            # Already in correct format
            pass
            
        data = {tag: swanlab.Image(img_tensor)}
        swanlab.log(data=data, step=int(global_step))

        return original_add_image(self, *args, **kwargs)
    
    def patched_add_text(self, *args, **kwargs):
        tag, text_string, global_step = _extract_args(
            args, kwargs, ['tag', 'text_string', 'global_step']
        )
        data = {tag: swanlab.Text(text_string)}
        swanlab.log(data=data, step=int(global_step))
        return original_add_text(self, *args, **kwargs)
    

    def patched_close(self):
        # 调用原始的close方法
        original_close(self)
        # 关闭SwanLab记录器
        swanlab.finish()

    return patched_init, patched_add_scalar, patched_add_scalars, patched_add_image, patched_add_text, patched_close


def _apply_patches(SummaryWriter, patched_methods):
    """
    应用monkey patch到SummaryWriter类
    
    Args:
        SummaryWriter: SummaryWriter类
        patched_methods: (patched_init, patched_add_scalar, patched_add_image, patched_close)
    """
    patched_init, patched_add_scalar, patched_add_scalars, patched_add_image, patched_add_text, patched_close = patched_methods
    
    SummaryWriter.__init__ = patched_init
    SummaryWriter.add_scalar = patched_add_scalar
    SummaryWriter.add_scalars = patched_add_scalars
    SummaryWriter.add_image = patched_add_image
    SummaryWriter.add_text = patched_add_text
    SummaryWriter.close = patched_close


def _sync_tensorboard_generic(import_func, logdir_extractor):
    """
    通用的tensorboard同步函数
    
    Args:
        import_func: 导入SummaryWriter的函数
        logdir_extractor: 提取logdir的函数
    """
    try:
        SummaryWriter = import_func()
    except ImportError as e:
        raise ImportError(f"Import failed: {e}")

    patched_methods = _create_patched_methods(SummaryWriter, logdir_extractor)
    _apply_patches(SummaryWriter, patched_methods)


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
    def import_tensorboardx():
        from tensorboardX import SummaryWriter
        return SummaryWriter
    
    def extract_logdir_tensorboardx(args, kwargs):
        logdir, _, _, _, _, _, _, log_dir, _ = _extract_args(
            args, kwargs, 
            ['logdir', 'comment', 'purge_step', 'max_queue', 'flush_secs', 
             'filename_suffix', 'write_to_disk', 'log_dir', 'comet_config']
        )
        return logdir or log_dir
    
    _sync_tensorboard_generic(import_tensorboardx, extract_logdir_tensorboardx)


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
    def import_torch_tensorboard():
        from torch.utils.tensorboard import SummaryWriter
        return SummaryWriter
    
    def extract_logdir_torch(args, kwargs):
        logdir, _ = _extract_args(args, kwargs, ['log_dir', 'comment'])
        return logdir
    
    _sync_tensorboard_generic(import_torch_tensorboard, extract_logdir_torch)
