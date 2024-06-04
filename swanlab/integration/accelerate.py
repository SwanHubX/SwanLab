"""
Docs: <WAITIGNG UPDATE>

For adaptation to the huggingface accelerate. You can used SwanLab as your tracker, experiment logs can be uploaded to
SwanLab or viewed using the local version of SwanLab. Detailed of used swanlab in accelerate train scripts are as follows:
------train.py in accelerate------
...
from swanlab.integration.accelerate import SwanLabTracker
...
tracker = SwanLabTracker("some_run_name")
accelerator = Accelerator(log_with=tracker)
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
            The name of the experiment run
        logging_dir (`str`, `os.PathLike`):
            Location for swanlab logs to be stored.
        kwargs:
            Additional key word arguments passed along to the `swanlab.init` method.
    """

    name = "swanlab"
    requires_logging_directory = False

    @on_main_process
    def __init__(self, run_name: str, logging_dir: Union[str, os.PathLike] = None, **kwargs):
        super().__init__()
        self.run_name = run_name
        self.logging_dir = os.path.join(logging_dir, run_name) if logging_dir is not None else None
        self.writer = swanlab.init(project=run_name, logdir=self.logging_dir, **kwargs)
        logger.debug(f"Initialized swanlab project {self.run_name} logging to {self.logging_dir}")
        logger.debug(
            "Make sure to log any initial configurations with `self.store_init_configuration` before training!"
        )

    @property
    def tracker(self):
        return self.writer

    @on_main_process
    def store_init_configuration(self, values: dict):
        """
        Logs `values` as hyperparameters for the run. Should be run at the beginning of your experiment. Stores the
        hyperparameters in a yaml file for future use.

        Args:
            values (Dictionary `str` to `bool`, `str`, `float` or `int`):
                Values to be stored as initial hyperparameters as key-value pairs. The values need to have type `bool`,
                `str`, `float`, `int`, or `None`.
        """
        swanlab.config.update(values)
        logger.debug("Stored initial configuration hyperparameters to swanlab")

    @on_main_process
    def log(self, values: dict, step: Optional[int] = None, **kwargs):
        """
        Logs `values` to the current run.

        Args:
            values (Dictionary `str` to `str`, `float`, `int` or `dict` of `str` to `float`/`int`):
                Values to be logged as key-value pairs. The values need to have type `str`, `float`, `int` or `dict` of
                `str` to `float`/`int`.
            step (`int`, *optional*):
                The run step. If included, the log will be affiliated with this step.
            kwargs:
                Additional key word arguments passed along to either `SummaryWriter.add_scaler`,
                `SummaryWriter.add_text`, or `SummaryWriter.add_scalers` method based on the contents of `values`.
        """
        self.writer.log(values, step=step, **kwargs)
        logger.debug("Successfully logged to swanlab")

    @on_main_process
    def log_images(self, values: dict, step: Optional[int], **kwargs):
        """
        Logs `images` to the current run.

        Args:
            values (Dictionary `str` to `List` of `np.ndarray` or `PIL.Image`):
                Values to be logged as key-value pairs. The values need to have type `List` of `np.ndarray` or
            step (`int`, *optional*):
                The run step. If included, the log will be affiliated with this step.
            kwargs:
                Additional key word arguments passed along to the `SummaryWriter.add_image` method.
        """
        for k, v in values.items():
            self.writer.add_images(k, v, global_step=step, **kwargs)
        logger.debug("Successfully logged images to swanlab")

    # @on_main_process
    # def log_table(
    #     self,
    #     table_name: str,
    #     columns: List[str] = None,
    #     data: List[List[Any]] = None,
    #     dataframe: Any = None,
    #     step: Optional[int] = None,
    #     **kwargs,
    # ):
    #     """
    #     Log a Table containing any object type (text, image, audio, video, molecule, html, etc). Can be defined either
    #     with `columns` and `data` or with `dataframe`.

    #     Args:
    #         table_name (`str`):
    #             The name to give to the logged table on the swanlab workspace
    #         columns (List of `str`'s *optional*):
    #             The name of the columns on the table
    #         data (List of List of Any data type *optional*):
    #             The data to be logged in the table
    #         dataframe (Any data type *optional*):
    #             The data to be logged in the table
    #         step (`int`, *optional*):
    #             The run step. If included, the log will be affiliated with this step.
    #     """

    #     values = {table_name: swanlab.Table(columns=columns, data=data, dataframe=dataframe)}
    #     self.log(values, step=step, **kwargs)

    @on_main_process
    def finish(self):
        """
        Closes `swanlab` writer
        """
        self.writer.finish()
        logger.debug("swanlab run closed")
