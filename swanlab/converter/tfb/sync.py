import functools
from typing import Optional, Sequence

import swanlab
from swanlab import vendor
from swanlab.converter.helper import extract_args
from swanlab.sdk.internal.pkg import safe


@safe.decorator(message="Failed to convert tensorboard image")
def _convert_tb_image(img_tensor, dataformats: str = "CHW"):
    """Convert a tensor/numpy image to HWC numpy array for swanlab.Image.

    Supports CHW, NCHW, HW, HWC dataformats.
    """
    np = vendor.np

    if hasattr(img_tensor, "cpu"):
        img_tensor = img_tensor.cpu()
    if hasattr(img_tensor, "numpy"):
        img_tensor = img_tensor.numpy()

    if dataformats == "CHW":
        img_tensor = np.transpose(img_tensor, (1, 2, 0))
    elif dataformats == "NCHW":
        img_tensor = np.transpose(img_tensor[0], (1, 2, 0))
    elif dataformats == "HW":
        img_tensor = np.expand_dims(img_tensor, axis=-1)

    return img_tensor


def _create_patched_methods(SummaryWriter, logdir_extractor, types=None):
    types_set = set(types) if types is not None else None

    original_init = SummaryWriter.__init__
    original_add_scalar = SummaryWriter.add_scalar
    original_add_scalars = SummaryWriter.add_scalars
    original_add_image = SummaryWriter.add_image
    original_add_text = SummaryWriter.add_text
    original_close = SummaryWriter.close

    def patched_init(self, *args, **kwargs):
        tb_logdir = logdir_extractor(args, kwargs)
        tb_config = {"tensorboard_logdir": tb_logdir}

        if not swanlab.has_run():
            swanlab.init(config=tb_config)
        else:
            swanlab.config.update(tb_config)

        return original_init(self, *args, **kwargs)

    @functools.wraps(original_add_scalar)
    def patched_add_scalar(self, *args, **kwargs):
        if types_set is not None and "scalar" not in types_set:
            return original_add_scalar(self, *args, **kwargs)

        tag, scalar_value, global_step = extract_args(args, kwargs, ["tag", "scalar_value", "global_step"])

        data = {tag: scalar_value}
        swanlab.log(data=data, step=int(global_step) if global_step is not None else None)

        return original_add_scalar(self, *args, **kwargs)

    @functools.wraps(original_add_scalars)
    def patched_add_scalars(self, *args, **kwargs):
        if types_set is not None and "scalars" not in types_set:
            return original_add_scalars(self, *args, **kwargs)

        tag, scalar_value_dict, global_step = extract_args(args, kwargs, ["tag", "scalar_value_dict", "global_step"])
        for dict_tag, value in scalar_value_dict.items():
            data = {f"{tag}/{dict_tag}": value}
            swanlab.log(data=data, step=int(global_step) if global_step is not None else None)

        return original_add_scalars(self, *args, **kwargs)

    @functools.wraps(original_add_image)
    def patched_add_image(self, *args, **kwargs):
        if types_set is not None and "image" not in types_set:
            return original_add_image(self, *args, **kwargs)

        tag, img_tensor, global_step, dataformats = extract_args(
            args, kwargs, ["tag", "img_tensor", "global_step", "dataformats"]
        )
        dataformats = dataformats or "CHW"

        converted = _convert_tb_image(img_tensor, dataformats)
        if converted is not None:
            data = {tag: swanlab.Image(converted)}
            swanlab.log(data=data, step=int(global_step) if global_step is not None else None)

        return original_add_image(self, *args, **kwargs)

    @functools.wraps(original_add_text)
    def patched_add_text(self, *args, **kwargs):
        if types_set is not None and "text" not in types_set:
            return original_add_text(self, *args, **kwargs)

        tag, text_string, global_step = extract_args(args, kwargs, ["tag", "text_string", "global_step"])
        data = {tag: swanlab.Text(text_string)}
        swanlab.log(data=data, step=int(global_step) if global_step is not None else None)

        return original_add_text(self, *args, **kwargs)

    def patched_close(self):
        original_close(self)
        swanlab.finish()

    return (patched_init, patched_add_scalar, patched_add_scalars, patched_add_image, patched_add_text, patched_close)


def _apply_patches(SummaryWriter, patched_methods):
    (patched_init, patched_add_scalar, patched_add_scalars, patched_add_image, patched_add_text, patched_close) = (
        patched_methods
    )

    SummaryWriter.__init__ = patched_init
    SummaryWriter.add_scalar = patched_add_scalar
    SummaryWriter.add_scalars = patched_add_scalars
    SummaryWriter.add_image = patched_add_image
    SummaryWriter.add_text = patched_add_text
    SummaryWriter.close = patched_close


def sync_tensorboardX(types: Optional[Sequence[str]] = None):
    """Monkey-patch tensorboardX SummaryWriter to forward logs to SwanLab.

    Call before creating a SummaryWriter instance.

    Args:
        types: Data types to sync. Options: 'scalar', 'scalars', 'image', 'text'.
               None syncs all types.

    Example::

        import swanlab
        swanlab.sync_tensorboardX()

        from tensorboardX import SummaryWriter
        writer = SummaryWriter('runs/example')
        writer.add_scalar('loss', 0.5, 0)
        writer.close()
    """
    vendor.tensorboardX
    from tensorboardX import SummaryWriter

    def extract_logdir(args, kwargs):
        logdir, _, _, _, _, _, _, log_dir, _ = extract_args(
            args,
            kwargs,
            [
                "logdir",
                "comment",
                "purge_step",
                "max_queue",
                "flush_secs",
                "filename_suffix",
                "write_to_disk",
                "log_dir",
                "comet_config",
            ],
        )
        return logdir or log_dir

    patched_methods = _create_patched_methods(SummaryWriter, extract_logdir, types)
    _apply_patches(SummaryWriter, patched_methods)


def sync_tensorboard_torch(types: Optional[Sequence[str]] = None):
    """Monkey-patch torch.utils.tensorboard SummaryWriter to forward logs to SwanLab.

    Call before creating a SummaryWriter instance.

    Args:
        types: Data types to sync. Options: 'scalar', 'scalars', 'image', 'text'.
               None syncs all types.

    Example::

        import swanlab
        swanlab.sync_tensorboard_torch()

        from torch.utils.tensorboard import SummaryWriter
        writer = SummaryWriter('runs/example')
        writer.add_scalar('loss', 0.5, 0)
        writer.close()
    """
    vendor.torch
    from torch.utils.tensorboard import SummaryWriter

    def extract_logdir(args, kwargs):
        log_dir, _ = extract_args(args, kwargs, ["log_dir", "comment"])
        return log_dir

    patched_methods = _create_patched_methods(SummaryWriter, extract_logdir, types)
    _apply_patches(SummaryWriter, patched_methods)
