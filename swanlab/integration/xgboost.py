import json
from typing import cast
import xgboost as xgb  # type: ignore
from xgboost import Booster
import swanlab


class SwanLabCallback(xgb.callback.TrainingCallback):
    def __init__(self):
        # 如果没有注册过实验
        swanlab.config["FRAMEWORK"] = "xgboost"
        if swanlab.get_run() is None:
            raise RuntimeError("You must call swanlab.init() before SwanLabCallback(). 你必须在SwanLabCallback()之前，调用swanlab.init().")
        
    def before_training(self, model: Booster) -> Booster:
        """Run before training is finished."""
        # Update SwanLab config
        config = model.save_config()
        swanlab.config.update(json.loads(config))

        return model

    def after_training(self, model: Booster) -> Booster:
        """Run after training is finished."""

        # Log the best score and best iteration
        if model.attr("best_score") is not None:
            swanlab.log(
                {
                    "best_score": float(cast(str, model.attr("best_score"))),
                    "best_iteration": int(cast(str, model.attr("best_iteration"))),
                }
            )

        return model

    def after_iteration(self, model: Booster, epoch: int, evals_log: dict) -> bool:
        """Run after each iteration. Return True when training should stop."""
        # Log metrics
        for data, metric in evals_log.items():
            for metric_name, log in metric.items():
                    swanlab.log({f"{data}-{metric_name}": log[-1]})

        swanlab.log({"epoch": epoch})

        return False

