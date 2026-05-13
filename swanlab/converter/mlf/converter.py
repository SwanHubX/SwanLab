from typing import Optional

import swanlab
from swanlab import vendor
from swanlab.converter.base import BaseConverter

mlflow = vendor.mlflow


class MLFlowConverter(BaseConverter):
    """Convert MLflow runs to SwanLab via the MLflow Tracking API.

    Mapping: MLflow experiment → SwanLab project, MLflow run → SwanLab run.

    Example::

        mlf = MLFlowConverter(project="my-swanlab-proj")
        mlf.run(tracking_uri="http://localhost:5000", experiment="my-experiment")
        mlf.run(tracking_uri="http://localhost:5000", run_id="abc123...")
    """

    def run(  # type: ignore[override]
        self,
        tracking_uri: str,
        experiment: Optional[str] = None,
        run_id: Optional[str] = None,
    ) -> None:

        client = mlflow.MlflowClient(tracking_uri=tracking_uri)

        # Resolve target runs
        if run_id is not None:
            target_runs = self._resolve_single_run(client, run_id)
            if target_runs is None:
                return
        else:
            target_runs = self._resolve_all_runs(client, experiment)
            if target_runs is None:
                return

        converted, failed = 0, 0
        for ex_name, run in target_runs:
            ok = self._convert_run(client, ex_name, run)
            if ok:
                converted += 1
            else:
                failed += 1

        if failed:
            print(f"\nConverted {converted} run(s), {failed} run(s) failed.")
        else:
            print(f"\nAll {converted} MLflow run(s) converted successfully.")

    def _resolve_single_run(self, client, run_id: str):
        """Resolve a single MLflow run by ID."""
        try:
            run = client.get_run(run_id)
        except Exception as e:
            print(f'Error: could not find MLflow run with id "{run_id}": {e}')
            return None
        ex = client.get_experiment(run.info.experiment_id)
        return [(ex.name, run)]

    def _resolve_all_runs(self, client, experiment: Optional[str]):
        """Resolve all runs from specified experiments (or all experiments)."""
        if experiment is None:
            experiments = client.search_experiments()
        else:
            try:
                ex = client.get_experiment(experiment)
            except Exception:
                ex = client.get_experiment_by_name(experiment)
            if not ex:
                print(f'Error: could not find experiment with id or name "{experiment}"')
                return None
            experiments = (ex,)

        result = []
        for ex in list(experiments)[::-1]:
            runs = client.search_runs(ex.experiment_id)
            print(f'Experiment "{ex.name}" — {len(runs)} run(s) to convert.')
            for run in runs:
                result.append((ex.name, run))
        return result

    def _convert_run(self, client, project_name: str, run) -> bool:
        """Convert a single MLflow run to SwanLab. Returns True on success."""
        run_id = run.info.run_id
        mlflow_run_name = run.data.tags.get("mlflow.runName")
        description = run.data.tags.get("mlflow.note.content")
        display_name = mlflow_run_name or run_id[:8]

        try:
            swanlab_run = swanlab.init(
                project=self.project or project_name,
                workspace=self.workspace,
                mode=self.mode,  # type: ignore[arg-type]
                log_dir=self.log_dir,
                name=mlflow_run_name,
                description=description,
                reinit=True,
            )
        except RuntimeError as e:
            print(f"  ✗ Failed to init run: {display_name}")
            print(f"    Error: {e}")
            return False

        try:
            swanlab_run.config.update(
                {
                    "mlflow_run_id": run_id,
                    "mlflow_run_name": mlflow_run_name,
                    "mlflow_run_description": description,
                    "mlflow_run_params": run.data.params,
                    "mlflow_run_tags": {k: v for k, v in run.data.tags.items() if not k.startswith("mlflow.")},
                }
            )

            all_metric_histories = {}
            for key in run.data.metrics.keys():
                all_metric_histories[key] = client.get_metric_history(run_id, key)

            total_metrics = sum(len(h) for h in all_metric_histories.values())
            print(f"  Converting run: {display_name} ({total_metrics} metrics)")

            for key, history in all_metric_histories.items():
                for m in history:
                    swanlab_run.log({m.key: m.value}, step=m.step)

            swanlab_run.finish()
            print(f"  ✓ Finished converting run: {display_name}")
            return True
        except Exception as e:
            print(f"  ✗ Failed to convert run: {display_name}")
            print(f"    Error: {e}")
            try:
                swanlab_run.finish()
            except Exception:
                pass
            return False
