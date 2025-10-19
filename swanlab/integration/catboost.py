try:
    import catboost
except ImportError:
    raise RuntimeError(
        "This module requires `catboost` to be installed. "
        "Please install it with command: pip install catboost"
    )

from typing import Any, Dict, List, NamedTuple, Optional
import swanlab


class IterationInfo(NamedTuple):
    iteration: int
    metrics: Dict[str, Dict[str, List[float]]]


class SwanLabCallback:
    def __init__(self, params: Optional[Dict[str, Any]] = None):
        """
        Initializes the SwanLabCallback for CatBoost.

        :params: Optional dictionary of parameters to log.
        """
        if swanlab.get_run() is None:
            raise RuntimeError(
                "You must call swanlab.init() before using SwanLabCallback."
            )
        swanlab.config["FRAMEWORK"] = "catboost"  # type: ignore
        if params and isinstance(params, dict):
            swanlab.config.update(params)

    def update_config(self, config: Dict[str, Any]):
        swanlab.config.update(config)

    def after_iteration(self, info: IterationInfo) -> bool:
        """
        Called after each iteration. Logs metrics and handles first-time setup.
        :return: True/False to continue training.
        """
        # Log metrics with data_name_metric_name format based on CatBoost JSON structure
        if info.metrics:
            for stage_name, metrics in info.metrics.items():
                for metric_name, values in metrics.items():
                    if values:
                        swanlab.log({f"{stage_name}_{metric_name}": values[-1]})

        return True
