"""
https://docs.swanlab.cn/guide_cloud/integration/integration-xgboost.html
"""

import json
from typing import cast, Any, Dict
import xgboost as xgb  # type: ignore
from xgboost import Booster
import swanlab


class SwanLabCallback(xgb.callback.TrainingCallback):
    def __init__(
        self,
        log_feature_importance: bool = True,
        importance_type: str = "gain",
        ):
        self.log_feature_importance = log_feature_importance
        self.importance_type = importance_type
        # 如果没有注册过实验
        swanlab.config["FRAMEWORK"] = "xgboost"
        if swanlab.get_run() is None:
            raise RuntimeError("You must call swanlab.init() before SwanLabCallback(). 你必须在SwanLabCallback()之前，调用swanlab.init().")
    
    def update_config(self, config: Dict[str, Any]):
        swanlab.config.update(config)
        
    def before_training(self, model: Booster) -> Booster:
        """Run before training is finished."""
        # Update SwanLab config
        config = model.save_config()
        swanlab.config.update(json.loads(config))

        return model

    def after_training(self, model: Booster) -> Booster:
        """Run after training is finished."""

        if self.log_feature_importance:
            self._log_feature_importance(model)
        
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

    def _log_feature_importance(self, model: Booster) -> None:
        fi = model.get_score(importance_type=self.importance_type)
        x = list(fi.keys())
        y = list(fi.values())
        y = [round(i, 2) for i in y]  # 保留两位小数
        bar = swanlab.echarts.Bar()
        bar.add_xaxis(x)
        bar.add_yaxis("Importance", y)
        swanlab.log(
            {
                "Feature Importance": bar
            }
        )