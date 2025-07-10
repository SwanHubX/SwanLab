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


def sync_wandb(
    mode:str="cloud",
    wandb_run:bool=True,
    workspace:str=None,
    logdir:str=None,
    ):
    """
    sync wandb with swanlab, 暂时不支持log非标量类型
    
    - mode: "cloud", "local" or "disabled". https://docs.swanlab.cn/api/py-init.html
    - wandb_run: 如果此参数设置为False，则不会将数据上传到wandb，等同于设置wandb.init(mode="offline")。
    - workspace: swanlab的组织空间username
    - logdir: swanlab日志文件存储目录
    
    usecase:
    ```python
    import wandb
    import random
    import swanlab

    swanlab.sync_wandb()
    # swanlab.init(project="sync_wandb")

    wandb.init(
        project="test",
        config={"a": 1, "b": 2},
        name="test",
    )

    epochs = 10
    offset = random.random() / 5
    for epoch in range(2, epochs):
        acc = 1 - 2 ** -epoch - random.random() / epoch - offset
        loss = 2 ** -epoch + random.random() / epoch + offset

        wandb.log({"acc": acc, "loss": loss})
    ```
    """
    try:
        import wandb
        from wandb import sdk as wandb_sdk
    except ImportError:
        raise ImportError("please install wandb first, command: `pip install wandb`")
    
    original_init = wandb.init
    original_log = wandb_sdk.wandb_run.Run.log
    original_finish = wandb_sdk.wandb_run.Run.finish
    original_config_update = wandb_sdk.wandb_config.Config.update
    
    def patched_init(*args, **kwargs):
        entity, project, dir, id, name, notes, tags, config, config_exclude_keys, reinit = _extract_args(
            args, kwargs, ['entity', 'project', 'dir', 'id', 'name', 'notes', 'tags', 'config', 'config_exclude_keys', 'reinit']
        )
        
        if swanlab.data.get_run() is None:
            swanlab.init(
                project=project,
                workspace=workspace,
                experiment_name=name,
                description=notes,
                config=config,
                tags=tags,
                mode=mode,
                logdir=logdir,
                reinit=reinit,
                )
        else:
            swanlab.config.update(config)
        
        if wandb_run is False:
            kwargs["mode"] = "offline"
            return original_init(*args, **kwargs)
        else:
            return original_init(*args, **kwargs)

    def patched_config_update(self, *args, **kwargs):
        d, _ = _extract_args(args, kwargs, ['d', 'allow_val_change'])
        
        if d is not None:
            swanlab.config.update(d)
        return original_config_update(self, *args, **kwargs)

    def patched_log(self, *args, **kwargs):
        data, step, commit, sync = _extract_args(args, kwargs, ['data', 'step', 'commit', 'sync'])
        
        if data is None:
            return original_log(self, *args, **kwargs)
        
        # 处理数据，支持 wandb.Image
        processed_data = {}
        for key, value in data.items():
            if isinstance(value, (int, float, bool, str)):
                # 标量类型直接保留
                processed_data[key] = value
            elif hasattr(value, '__class__') and value.__class__.__name__ == 'Image' and hasattr(value, 'image'):
                # 检测是否为 wandb.Image
                try:
                    # 获取 wandb.Image 的图像数据
                    if value.image is not None:
                        # 将 PIL Image 转换为 numpy 数组
                        import numpy as np
                        img_array = np.array(value.image)
                        
                        # 创建 swanlab.Image
                        caption = getattr(value, '_caption', None)
                        swanlab_image = swanlab.Image(img_array, caption=caption)
                        processed_data[key] = swanlab_image
                    else:
                        # 如果 image 为 None，尝试使用 _image
                        if hasattr(value, '_image') and value._image is not None:
                            import numpy as np
                            img_array = np.array(value._image)
                            caption = getattr(value, '_caption', None)
                            swanlab_image = swanlab.Image(img_array, caption=caption)
                            processed_data[key] = swanlab_image
                except Exception as e:
                    # 如果转换失败，记录错误但继续处理其他数据
                    print(f"Warning: Failed to convert wandb.Image for key '{key}': {e}")
                    continue
            elif isinstance(value, list) and value and hasattr(value[0], '__class__') and value[0].__class__.__name__ == 'Image':
                # 检测是否为 wandb.Image 列表
                try:
                    import numpy as np
                    swanlab_images = []
                    for v in value:
                        if hasattr(v, 'image') and v.image is not None:
                            img_array = np.array(v.image)
                            caption = getattr(v, '_caption', None)
                            swanlab_images.append(swanlab.Image(img_array, caption=caption))
                        elif hasattr(v, '_image') and v._image is not None:
                            img_array = np.array(v._image)
                            caption = getattr(v, '_caption', None)
                            swanlab_images.append(swanlab.Image(img_array, caption=caption))
                    if swanlab_images:
                        processed_data[key] = swanlab_images
                except Exception as e:
                    # 如果转换失败，记录错误但继续处理其他数据
                    print(f"Warning: Failed to convert wandb.Image list for key '{key}': {e}")
                    continue
        
        if processed_data:
            swanlab.log(data=processed_data, step=step)
        
        return original_log(self, *args, **kwargs)
    
    def patched_finish(*args, **kwargs):
        swanlab.finish()
        return original_finish(*args, **kwargs)

    wandb.init = patched_init
    wandb_sdk.wandb_run.Run.log = patched_log
    wandb_sdk.wandb_run.Run.finish = patched_finish
    wandb_sdk.wandb_config.Config.update = patched_config_update