from typing import Optional

from swanlab import vendor
from swanlab.converter.base import BaseConverter


class MLFlowConverter(BaseConverter):
    """Convert MLflow runs to SwanLab via the MLflow Tracking API.

    Example::

        mlf = MLFlowConverter(project="my-swanlab-proj")
        mlf.run(tracking_uri="http://localhost:5000", experiment="my-experiment")
    """

    def run(  # type: ignore[override]
        self,
        tracking_uri: str,
        experiment: Optional[str] = None,
    ) -> None:
        mlflow = vendor.mlflow
        import swanlab

        client = mlflow.MlflowClient(tracking_uri=tracking_uri)

        # Resolve experiments
        if experiment is None:
            experiments = client.search_experiments()
        else:
            try:
                ex = client.get_experiment(experiment)
            except mlflow.exceptions.MlflowException:
                ex = client.get_experiment_by_name(experiment)
            if not ex:
                print(f'Error: could not find experiment with id or name "{experiment}"')
                return
            experiments = (ex,)

        # Reverse so oldest experiment is converted first
        experiments = list(experiments)[::-1]

        for ex in experiments:
            runs = client.search_runs(ex.experiment_id)
            print(f'Experiment "{ex.name}" — {len(runs)} run(s) to convert.')

            for run in runs:
                run_id = run.info.run_id
                mlflow_run_name = run.data.tags.get("mlflow.runName")
                description = run.data.tags.get("mlflow.note.content")

                swanlab_run = swanlab.init(
                    project=self.project or ex.name,
                    workspace=self.workspace,
                    mode=self.mode,  # type: ignore[arg-type]
                    log_dir=self.log_dir,
                    name=mlflow_run_name,
                    description=description,
                    reinit=True,
                )

                swanlab_run.config.update(
                    {
                        "mlflow_run_id": run_id,
                        "mlflow_run_name": mlflow_run_name,
                        "mlflow_run_description": description,
                        "mlflow_run_params": run.data.params,
                        "mlflow_run_tags": {k: v for k, v in run.data.tags.items() if not k.startswith("mlflow.")},
                    }
                )

                # Prefetch all metric histories to avoid repeated API calls
                all_metric_histories = {}
                for key in run.data.metrics.keys():
                    all_metric_histories[key] = client.get_metric_history(run_id, key)

                total_metrics = sum(len(h) for h in all_metric_histories.values())
                print(f"  Converting run: {mlflow_run_name or run_id[:8]} ({total_metrics} metrics)")

                for key, history in all_metric_histories.items():
                    for m in history:
                        swanlab_run.log({m.key: m.value}, step=m.step)

                swanlab_run.finish()
                print(f"  ✓ Finished converting run: {mlflow_run_name or run_id[:8]}")

        print("All MLflow runs converted successfully.")
