"""
Docs: https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-accelerate.html

For adaptation to the huggingface accelerate. You can used SwanLab as your tracker, experiment logs can be uploaded to
SwanLab or viewed using the local version of SwanLab. Detailed of used swanlab in accelerate train scripts are as follows:
------train.py in accelerate------
...
from swanlab.integration.accelerate import SwanLabTracker
...
tracker = SwanLabTracker("some_run_name")
accelerator = Accelerator(log_with=[tracker])
...
---------------------------------
These also can be mixed with existing trackers, including with "all":
------train.py in accelerate------
...
from swanlab.integration.accelerate import SwanLabTracker
...
tracker = SwanLabTracker("some_run_name")
accelerator = Accelerator(log_with=[tracker, "all"])
...
---------------------------------
"""

import os
from typing import Any, Dict, List, Optional, Union

import swanlab

try:
    from accelerate.tracking import GeneralTracker, on_main_process
    from accelerate.logging import get_logger
except ImportError:
    raise RuntimeError(
        "This contrib module requires Accelerate to be installed. "
        "Please install it with command: \n pip install accelerate"
    )


logger = get_logger(__name__)


class SwanLabTracker(GeneralTracker):
    """
    A `Tracker` class that supports `swanlab`. Should be initialized at the start of your script.

    Args:
        run_name (`str`):
            The name of the experiment run.
        **kwargs (additional keyword arguments, *optional*):
            Additional key word arguments passed along to the `swanlab.init` method.
    """

    name = "swanlab"
    requires_logging_directory = False
    main_process_only = True  # it must be True in tracker module

    def __init__(self, run_name: str, **kwargs):
        super().__init__()
        self.run_name = run_name
        self.init_kwargs = kwargs
        self.start()  # auto start swanlab when instantiation tracker

    @on_main_process
    def start(self):
        if hasattr(self, "run") and self.run is not None:
            return

        self.run = swanlab.init(project=self.run_name, **self.init_kwargs)
        swanlab.config["FRAMEWORK"] = "accelerate"  # add accelerate logo in config
        logger.debug(f"Initialized SwanLab project {self.run_name}")
        logger.debug(
            "Make sure to log any initial configurations with `self.store_init_configuration` before training!"
        )

    @property
    def tracker(self):
        return self.run

    @on_main_process
    def store_init_configuration(self, values: dict):
        """
        Logs `values` as hyperparameters for the run. Should be run at the beginning of your experiment.

        Args:
            values (Dictionary `str` to `bool`, `str`, `float` or `int`):
                Values to be stored as initial hyperparameters as key-value pairs. The values need to have type `bool`,
                `str`, `float`, `int`, or `None`.
        """
        import swanlab

        swanlab.config.update(values, allow_val_change=True)
        logger.debug("Stored initial configuration hyperparameters to SwanLab")

    @on_main_process
    def log(self, values: dict, step: Optional[int] = None, **kwargs):
        """
        Logs `values` to the current run.

        Args:
        data : Dict[str, DataType]
            Data must be a dict.
            The key must be a string with 0-9, a-z, A-Z, " ", "_", "-", "/".
            The value must be a `float`, `float convertible object`, `int` or `swanlab.data.BaseType`.
        step : int, optional
            The step number of the current data, if not provided, it will be automatically incremented.
        If step is duplicated, the data will be ignored.
            kwargs:
                Additional key word arguments passed along to the `swanlab.log` method. Likes:
                    print_to_console : bool, optional
                        Whether to print the data to the console, the default is False.
        """
        self.run.log(values, step=step, **kwargs)
        logger.debug("Successfully logged to SwanLab")

    @on_main_process
    def log_images(self, values: dict, step: Optional[int] = None, **kwargs):
        """
        Logs `images` to the current run.

        Args:
            values (Dictionary `str` to `List` of `np.ndarray` or `PIL.Image`):
                Values to be logged as key-value pairs. The values need to have type `List` of `np.ndarray` or
            step (`int`, *optional*):
                The run step. If included, the log will be affiliated with this step.
            kwargs:
                Additional key word arguments passed along to the `swanlab.log` method. Likes:
                    print_to_console : bool, optional
                        Whether to print the data to the console, the default is False.
        """
        import swanlab

        for k, v in values.items():
            self.log({k: [swanlab.Image(image) for image in v]}, step=step, **kwargs)
        logger.debug("Successfully logged images to SwanLab")

    @on_main_process
    def finish(self):
        """
        Closes `swanlab` writer
        """
        self.run.finish()
        logger.debug("SwanLab run closed")
