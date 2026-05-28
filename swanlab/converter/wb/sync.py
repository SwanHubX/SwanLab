from typing import Optional

import swanlab
from swanlab import vendor
from swanlab.converter.helper import extract_args
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.run import ModeType


def sync_wandb(
    mode: ModeType = "online",
    wandb_run: bool = True,
    workspace: Optional[str] = None,
    log_dir: Optional[str] = None,
):
    """Monkey-patch wandb，将日志实时转发到 SwanLab。

    在调用 ``wandb.init()`` 之前调用此函数。

    Args:
        mode: SwanLab mode — online | local | offline | disabled
        wandb_run: False 时将 wandb 设为 offline（不实际上传 wandb）
        workspace: SwanLab workspace
        log_dir: SwanLab 日志目录

    Example::

        import swanlab
        swanlab.sync_wandb()

        import wandb
        wandb.init(project="test", config={"lr": 0.01})
        wandb.log({"loss": 0.5})
        wandb.finish()
    """
    wandb = vendor.wandb
    WandbImage = wandb.Image
    wandb_sdk = wandb.sdk  # type: ignore

    original_init = wandb.init
    original_log = wandb_sdk.wandb_run.Run.log
    original_finish = wandb_sdk.wandb_run.Run.finish
    original_config_update = wandb_sdk.wandb_config.Config.update

    def patched_init(*args, **kwargs):
        (
            _entity,
            project,
            _dir,
            _run_id,
            name,
            notes,
            tags,
            config,
            _config_exclude_keys,
            reinit,
            group,
            job_type,
        ) = extract_args(
            original_init,
            args,
            kwargs,
            [
                "entity",
                "project",
                "dir",
                "id",
                "name",
                "notes",
                "tags",
                "config",
                "config_exclude_keys",
                "reinit",
                "group",
                "job_type",
            ],
        )

        if not swanlab.has_run():
            swanlab.init(
                project=project,
                workspace=workspace,
                name=name,
                description=notes,
                config=config,
                tags=tags,
                mode=mode,
                log_dir=log_dir,
                reinit=reinit,
                group=group,
                job_type=job_type,
            )
        else:
            swanlab.config.update(config or {})

        if wandb_run is False:
            kwargs["mode"] = "offline"

        return original_init(*args, **kwargs)

    def patched_config_update(self, *args, **kwargs):
        d, _ = extract_args(original_config_update, (self,) + args, kwargs, ["d", "allow_val_change"])
        if d is not None:
            swanlab.config.update(d)

        extra_kwargs = {k: v for k, v in kwargs.items() if k not in {"d", "allow_val_change"}}
        if extra_kwargs:
            swanlab.config.update(extra_kwargs)

        return original_config_update(self, *args, **kwargs)

    def patched_log(self, *args, **kwargs):
        data, step, _commit, _sync = extract_args(
            original_log,
            (self,) + args,
            kwargs,
            ["data", "step", "commit", "sync"],
        )

        if data is None or not hasattr(data, "items"):
            return original_log(self, *args, **kwargs)

        processed_data = {}
        for key, value in data.items():
            if isinstance(value, (int, float, bool, str)):
                processed_data[key] = value
            elif isinstance(value, WandbImage):
                processed_data[key] = _convert_wandb_image(value)
            elif isinstance(value, list) and value and isinstance(value[0], WandbImage):
                images = [_convert_wandb_image(v) for v in value]
                images = [img for img in images if img is not None]
                if images:
                    processed_data[key] = images

        processed_data = {k: v for k, v in processed_data.items() if v is not None}

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


@safe.decorator(message="Failed to convert wandb.Image")
def _convert_wandb_image(wandb_img):
    """将单个 wandb.Image 转换为 swanlab.Image，失败返回 None。"""
    np = vendor.np
    pil_image = getattr(wandb_img, "image", None) or getattr(wandb_img, "_image", None)
    if pil_image is None:
        return None
    img_array = np.array(pil_image)
    caption = getattr(wandb_img, "_caption", None)
    return swanlab.Image(img_array, caption=caption)
