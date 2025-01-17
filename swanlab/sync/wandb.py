import swanlab

def sync_wandb():
    """
    sync wandb with swanlab
    
    usecase:
    ```python
    import swanlab
    swanlab.sync_wandb()
    
    wandb.init(
        project="test",
        config={"a": 1, "b": 2},
        name="test",
    )

    for epoch in range(10):
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
    
    def patched_init(*args, **kwargs):
        project = kwargs.get('project', None)
        name = kwargs.get('name', None)
        config = kwargs.get('config', None)
        
        if swanlab.data.get_run() is None:
            swanlab.init(
                project=project,
                experiment_name=name,
                config=config)
        else:
            swanlab.config.update(config)
        
        return original_init(*args, **kwargs)

    def patched_log(*args, **kwargs):
        data = args[1]
        step = kwargs.get('step', None)
        
        # 过滤掉非标量类型
        filtered_data = {}
        for key, value in data.items():
            if isinstance(value, (int, float, bool, str)):
                filtered_data[key] = value
        
        swanlab.log(data=filtered_data, step=step)
        
        return original_log(*args, **kwargs)
    
    def patched_finish(*args, **kwargs):
        swanlab.finish()
        return original_finish(*args, **kwargs)

    wandb.init = patched_init
    wandb_sdk.wandb_run.Run.log = patched_log
    wandb_sdk.wandb_run.Run.finish = patched_finish