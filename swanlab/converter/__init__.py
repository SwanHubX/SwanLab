"""
Public converter API — convert experiment logs from other tracking tools into SwanLab.

Usage::

    from swanlab.converter import WandbConverter

    wb = WandbConverter(
        project="my-project",   # optional SwanLab project name
        mode="online",          # online | local | offline | disabled
    )
    wb.run(wb_project="WANDB_PROJECT", wb_entity="WANDB_ENTITY")

For local wandb files::

    from swanlab.converter import WandbLocalConverter

    wb = WandbLocalConverter(project="my-project")
    wb.run(root_wandb_dir="./wandb")

For TensorBoard::

    from swanlab.converter import TFBConverter

    tb_converter = TFBConverter(convert_dir="[TFEVENT_LOGDIR]")
    tb_converter.run()

For MLflow::

    from swanlab.converter import MLFlowConverter

    mlf = MLFlowConverter(project="my-project")
    # Optional: `experiment`
    mlf.run(tracking_uri="http://localhost:5000", experiment: Optional[str]="1")
"""

from swanlab.converter.mlf import MLFlowConverter
from swanlab.converter.tfb import TFBConverter
from swanlab.converter.wb import WandbConverter, WandbLocalConverter

__all__ = ["WandbConverter", "WandbLocalConverter", "TFBConverter", "MLFlowConverter"]
