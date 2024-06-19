"""
Docs: https://docs.swanlab.cn/zh/guide_cloud/integration/integration-mmengine.html

For adaptation to the mmengine framework, this adaptation also applies to frameworks such as mmdetection, xtuner, etc.
, which use mmengine as the engine. By setting 'vis_backends' to 'swanlab' in the config file, experiment logs can be
uploaded to SwanLab or viewed using the local version of SwanLab. Detailed configuration file changes are as follows:
------config.py in mmengine------
...
custom_imports = dict(
    imports=["swanlab.integration.mmengine"], allow_failed_imports=False
)
vis_backends = [
    dict(
        type="SwanlabVisBackend",
        save_dir="./swanlog_save_path",             # swanlab save path
        init_kwargs={                               # swanlab.init args
            "project": "YourProject ",              # project name on swanlab
            "experiment_name": "YourExperiment",    # experiment name on swanlab
            "description": "have fun",              # experiment description (can be null)
            "workspace": "YourOrganization",        # Your Organization on swanlab
            # "mode": "cloud",                       # Upload to cloud
        },
    ),
]
...
---------------------------------
"""

import os
from typing import Any, Callable, List, Optional, Sequence, Union

import numpy as np
import torch

try:
    import mmengine
except:
    raise ValueError(
        "This module requires mmengine to be installed. Please install it with command: \n pip install mmengine."
        "\nMore details can be found at https://mmengine.readthedocs.io/en/latest/get_started/installation.html"
    )

from mmengine.registry import VISBACKENDS
from mmengine.visualization.vis_backend import BaseVisBackend, force_init_env
from mmengine.config import Config


@VISBACKENDS.register_module()
class SwanlabVisBackend(BaseVisBackend):
    """Swanlab visualization backend class for mmengine.
    Examples:
        >>> from mmengine.visualization import SwanlabVisBackend
        >>> import numpy as np
        >>> Swanlab_vis_backend = SwanlabVisBackend()
        >>> img=np.random.randint(0, 256, size=(10, 10, 3))
        >>> Swanlab_vis_backend.add_image('img', img)
        >>> Swanlab_vis_backend.add_scaler('mAP', 0.6)
        >>> Swanlab_vis_backend.add_scalars({'loss': [1, 2, 3],'acc': 0.8})
        >>> cfg = Config(dict(a=1, b=dict(b1=[0, 1])))
        >>> Swanlab_vis_backend.add_config(cfg)

    Args:
        save_dir (str, optional): The root directory to save the files
            produced by the visualizer. Default used './swanlab'
        init_kwargs (dict, optional): Swanlab initialization
            input parameters.
            See `swanlab.init <NEED UPDATE>`_ for
            details. Defaults to None.
    """

    def __init__(
        self,
        save_dir: str = None,
        init_kwargs: Optional[dict] = None,
    ):
        self._save_dir = save_dir
        self._env_initialized = False
        self._init_kwargs = init_kwargs

    @force_init_env
    def experiment(self) -> Any:
        """Return the experiment object associated with this visualization
        backend.

        The experiment attribute can get the swanlab backend. If you want
        to write other data, such as writing a table, you can directly get
        the visualization backend through experiment.
        """
        return self._swanlab

    def _init_env(self) -> Any:
        """Setup env for swanlab."""
        if self._save_dir is not None:
            if not os.path.exists(self._save_dir):
                os.makedirs(self._save_dir, exist_ok=True)  # type: ignore
        if self._init_kwargs is None:
            self._init_kwargs = {"logdir": self._save_dir}
        else:
            self._init_kwargs.setdefault("logdir", self._save_dir)
        try:
            import swanlab
        except ImportError:
            raise ImportError('Please run "pip install swanlab" to install swanlab')

        swanlab.init(**self._init_kwargs)
        self._swanlab = swanlab

    @force_init_env
    def add_config(self, config: Config, **kwargs) -> None:
        """Record the config to swanlab.

        Args:
            config (Config): The Config object
        """

        def repack_dict(a, prefix=""):
            """
            Unpack Nested Dictionary func
            """
            new_dict = dict()
            for key, value in a.items():
                key = str(key)
                if isinstance(value, dict):
                    if prefix != "":
                        new_dict.update(repack_dict(value, f"{prefix}/{key}"))
                    else:
                        new_dict.update(repack_dict(value, key))
                elif isinstance(value, list) or isinstance(value, tuple):
                    if all(not isinstance(element, dict) for element in value):
                        new_dict[key] = value
                    else:
                        for i, item in enumerate(value):
                            new_dict.update(repack_dict(item, f"{key}[{i}]"))
                elif prefix != "":
                    new_dict[f"{prefix}/{key}"] = value
                else:
                    new_dict[key] = value
            return new_dict

        config_dict = config.to_dict()
        self._swanlab.config.update(repack_dict(config_dict))

    @force_init_env
    def add_graph(self, model: torch.nn.Module, data_batch: Sequence[dict], **kwargs) -> None:
        """Record the model graph to swanlab.

        Args:
            model (torch.nn.Module): Model to draw.
            data_batch (Sequence[dict]): Batch of data from dataloader.
        """
        # todo: waiting for update
        pass

    @force_init_env
    def add_image(self, name: str, image: np.ndarray, step: int = 0, **kwargs) -> None:
        """Record the image to swanlab.

        Args:
            name (str): The image identifier.
            image (np.ndarray): The image to be saved. The format
                should be RGB. Defaults to None.
            step (int): Global step value to record. Defaults to 0.
        """
        image = self._swanlab.Image(image)
        self._swanlab.log({name: image}, step=step)

    @force_init_env
    def add_scalar(self, name: str, value: Union[int, float], step: int = 0, **kwargs) -> None:
        """Record the scalar to swanlab.

        Args:
            name (str): The scalar identifier.
            value (int, float): Value to save.
            step (int): Global step value to record. Defaults to 0.
        """
        self._swanlab.log({name: value}, step=step)

    @force_init_env
    def add_scalars(
        self,
        scalar_dict: dict,
        step: int = 0,
        file_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        """Record the scalars' data.

        Args:
            scalar_dict (dict): Key-value pair storing the tag and
                corresponding values.
            step (int): Global step value to record. Defaults to 0.
            file_path (str, optional): The scalar's data will be
                saved to the `file_path` file at the same time
                if the `file_path` parameter is specified.
                Defaults to None.
        """
        self._swanlab.log(scalar_dict, step=step)

    def close(self) -> None:
        """close an opened swanlab object."""
        if hasattr(self, "_swanlab"):
            self._swanlab.finish()
