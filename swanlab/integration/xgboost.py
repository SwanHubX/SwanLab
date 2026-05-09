from __future__ import annotations

import json
from typing import TYPE_CHECKING, cast

import swanlab
import swanlab.vendor
from swanlab import Callback

if TYPE_CHECKING:
    from typing import Any


class SwanLabCallback(Callback, swanlab.vendor.xgboost.callback.TrainingCallback):
    def __init__(
        self,
        log_feature_importance: bool = True,
        importance_type: str = "gain",
    ):
        self._log_feature_importance = log_feature_importance
        self._importance_type = importance_type

    @property
    def name(self) -> str:
        return "swanlab-integration-xgboost"

    def before_training(self, model: Any) -> Any:
        run = swanlab.get_run()
        if run is None:
            return model
        run.config["FRAMEWORK"] = "xgboost"
        config = model.save_config()
        run.config.update(json.loads(config))
        return model

    def after_training(self, model: Any) -> Any:
        if self._log_feature_importance:
            self._log_feature_importance_chart(model)

        if model.attr("best_score") is not None:
            swanlab.log(
                {
                    "best_score": float(cast(str, model.attr("best_score"))),
                    "best_iteration": int(cast(str, model.attr("best_iteration"))),
                },
            )

        return model

    def after_iteration(self, model: Any, epoch: int, evals_log: dict) -> bool:
        for data, metric in evals_log.items():
            for metric_name, log in metric.items():
                swanlab.log({f"{data}-{metric_name}": log[-1]}, step=epoch)

        return False

    def _log_feature_importance_chart(self, model: Any) -> None:
        fi = model.get_score(importance_type=self._importance_type)
        if not fi:
            return
        x = list(fi.keys())
        y = [round(v, 2) for v in fi.values()]
        bar = swanlab.echarts.Bar()
        bar.add_xaxis(x)
        bar.add_yaxis("Importance", y)
        swanlab.log({"Feature Importance": bar})
