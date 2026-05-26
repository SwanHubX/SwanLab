import os

import swanlab
from swanlab import vendor
from swanlab.converter.helper import extract_args
from swanlab.sdk.typings.run import ModeType


def sync_mlflow(mode: ModeType = "online") -> None:
    """Monkey-patch mlflow to forward logs to SwanLab.

    Call before any mlflow API usage.

    Args:
        mode: SwanLab mode — online | local | offline | disabled.

    Example::

        import swanlab
        swanlab.sync_mlflow()

        import mlflow
        mlflow.set_experiment("my_experiment")

        with mlflow.start_run(run_name="my_run"):
            mlflow.log_param("learning_rate", 0.01)
            for step in range(10):
                mlflow.log_metric("loss", 0.5 / (step + 1), step=step)

        mlflow.end_run()
    """
    mlflow = vendor.mlflow

    original_set_experiment = mlflow.set_experiment
    original_start_run = mlflow.start_run
    original_end_run = mlflow.end_run
    original_log_param = mlflow.log_param
    original_log_params = mlflow.log_params
    original_log_metric = mlflow.log_metric
    original_log_metrics = mlflow.log_metrics

    def patched_set_experiment(experiment_name, *args, **kwargs):
        if not swanlab.has_run():
            os.environ["SWANLAB_MLFLOW_PROJECT"] = experiment_name
        return original_set_experiment(experiment_name, *args, **kwargs)

    def patched_start_run(*args, **kwargs):
        run_name = extract_args(args, kwargs, ["run_name"])[0]

        if not swanlab.has_run():
            project = os.environ.get("SWANLAB_MLFLOW_PROJECT")
            swanlab.init(project=project, name=run_name, mode=mode)
        elif run_name:
            swanlab.config.update({"mlflow_run_name": run_name})

        return original_start_run(*args, **kwargs)

    def patched_end_run(*args, **kwargs):
        swanlab.finish()
        return original_end_run(*args, **kwargs)

    def patched_log_param(key, value, *args, **kwargs):
        swanlab.config.update({key: value})
        return original_log_param(key, value, *args, **kwargs)

    def patched_log_params(params, *args, **kwargs):
        swanlab.config.update(params)
        return original_log_params(params, *args, **kwargs)

    def patched_log_metric(key, value, *args, **kwargs):
        step = extract_args(args, kwargs, ["step"])[0]
        swanlab.log({key: value}, step=step)
        return original_log_metric(key, value, *args, **kwargs)

    def patched_log_metrics(metrics, *args, **kwargs):
        step = extract_args(args, kwargs, ["step"])[0]
        swanlab.log(metrics, step=step)
        return original_log_metrics(metrics, *args, **kwargs)

    mlflow.set_experiment = patched_set_experiment
    mlflow.start_run = patched_start_run
    mlflow.end_run = patched_end_run
    mlflow.log_param = patched_log_param
    mlflow.log_params = patched_log_params
    mlflow.log_metric = patched_log_metric
    mlflow.log_metrics = patched_log_metrics
