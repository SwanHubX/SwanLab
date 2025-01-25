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


def sync_wandb(mode:str="cloud", wandb_run:bool=True):
    """
    sync wandb with swanlab, 暂时不支持log非标量类型
    
    - mode: "cloud", "local" or "disabled". https://docs.swanlab.cn/api/py-init.html
    - wandb_run: 如果此参数设置为False，则不会将数据上传到wandb，等同于设置wandb.init(mode="offline")。
    
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
    original_finish = wandb_sdk.finish
    original_config_update = wandb_sdk.wandb_config.Config.update
    
    def patched_init(*args, **kwargs):
        entity, project, dir, id, name, notes, tags, config, config_exclude_keys = _extract_args(
            args, kwargs, ['entity', 'project', 'dir', 'id', 'name', 'notes', 'tags', 'config', 'config_exclude_keys']
        )
        
        if swanlab.data.get_run() is None:
            swanlab.init(
                project=project,
                experiment_name=name,
                description=notes,
                config=config,
                mode=mode)
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
        
        # 过滤掉非标量类型
        filtered_data = {}
        for key, value in data.items():
            if isinstance(value, (int, float, bool, str)):
                filtered_data[key] = value
        
        swanlab.log(data=filtered_data, step=step)
        
        return original_log(self, *args, **kwargs)
    
    def patched_finish(*args, **kwargs):
        swanlab.finish()
        return original_finish(*args, **kwargs)

    wandb.init = patched_init
    wandb_sdk.wandb_run.Run.log = patched_log
    wandb_sdk.wandb_run.Run.finish = patched_finish
    wandb_sdk.wandb_config.Config.update = patched_config_update