try:
    import catboost
except ImportError:
    raise RuntimeError(
        "This module requires `catboost` to be installed. " "Please install it with command: pip install catboost"
    )

from typing import Any, Dict
import swanlab


class SwanLabCallback:
    def __init__(self, log_params: bool = True, log_feature_importance: bool = True):
        """
        Initializes the SwanLabCallback for CatBoost.

        :param log_params: If True, logs the model's parameters.
        :param log_feature_importance: If True, logs feature importance after training.
        """
        if swanlab.get_run() is None:
            raise RuntimeError("You must call swanlab.init() before using SwanLabCallback.")
        self.log_params = log_params
        self.log_feature_importance = log_feature_importance
        self.model = None
        self.is_first_iteration = True

    def update_config(self, config: Dict[str, Any]):
        """Update SwanLab config with additional parameters."""
        swanlab.config.update(config)

    def after_iteration(self, info: Any) -> bool:
        """
        Called after each iteration. Logs metrics and handles first-time setup.
        :return: True to continue training.
        """
        # On first iteration, log framework and parameters
        if self.is_first_iteration:
            swanlab.config["FRAMEWORK"] = "catboost"
            if self.log_params:
                try:
                    # Try to extract model parameters from the info object
                    if hasattr(info, 'params') and info.params:
                        swanlab.config.update(info.params)
                    # Note: CatBoost doesn't provide direct access to model params in callbacks
                    # Users should use update_config method to log custom parameters
                except Exception:
                    print("swanlab: Could not log parameters for catboost. Please report this issue.")
            self.is_first_iteration = False

        # Log metrics
        if info.metrics:
            for metric_name, values in info.metrics.items():
                for dataset_name, metric_values in values.items():
                    swanlab.log({f"{dataset_name}_{metric_name}": metric_values[-1]})
        swanlab.log({"iteration": info.iteration})
        
        # Store reference to model for feature importance logging
        if hasattr(info, 'model'):
            self.model = info.model
        
        return True

    def log_feature_importance(self, model=None):
        """
        Log feature importance. This method should be called manually after training.
        
        :param model: Trained CatBoost model. If None, uses the model from the last callback.
        """
        target_model = model if model is not None else self.model
        if target_model is None:
            print("swanlab: No model available for feature importance logging.")
            return
            
        try:
            # Get feature importance
            feature_importance = target_model.get_feature_importance()
            feature_names = target_model.feature_names_
            
            # Create bar chart
            if feature_importance is not None and feature_names is not None:
                # Round importance values for better readability
                importance_values = [round(float(val), 2) for val in feature_importance]
                
                # Create bar chart using SwanLab
                bar = swanlab.echarts.Bar()
                bar.add_xaxis(list(feature_names))
                bar.add_yaxis("Importance", importance_values)
                swanlab.log({"Feature Importance": bar})
        except Exception as e:
            print(f"swanlab: Error logging feature importance: {e}")

    def log_best_metrics(self, best_score=None, best_iteration=None):
        """
        Log best metrics. This method should be called manually after training.
        
        :param best_score: Best score achieved during training.
        :param best_iteration: Best iteration during training.
        """
        if best_score is not None:
            swanlab.log({"best_score": best_score})
        
        if best_iteration is not None:
            swanlab.log({"best_iteration": best_iteration})