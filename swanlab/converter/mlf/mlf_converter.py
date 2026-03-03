"""
------example.py------
from swanlab.converter import MLFlowConverter

mlf_converter = MLFlowConverter()
mlf_converter.run(tracking_uri="MLFLOW_TRACKING_URL", experiment=0)
"""
import swanlab
from swanlab.log import swanlog as swl
import time

try:
    from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
    from rich.progress import ProgressColumn
    from rich.text import Text
except ImportError:
    Progress = None
    ProgressColumn = None
    Text = None

class MLFLowConverter:
    def __init__(
        self,
        project: str = None,
        workspace: str = None,
        mode: str = "cloud",
        logdir: str = None,
        
        **kwargs,
    ):
        
        self.project = project
        self.workspace = workspace
        self.mode = mode
        self.logdir = logdir
    
    def parse_mlflow_logs(self, tracking_uri:str, experiment:int = None):
    
        try:
            import mlflow
        except ImportError as e:
            raise TypeError(
                "MLFlow Converter requires mlflow. Install with 'pip install mlflow'."
            )

        client = mlflow.MlflowClient(tracking_uri=tracking_uri)

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

        # Reverse the order of experiments
        if isinstance(experiments, tuple):
            experiments = experiments[::-1]
        else:
            experiments = list(experiments)[::-1]

        for ex in experiments:
            runs = client.search_runs(ex.experiment_id)
            for run in runs:
                run_id = run.info.run_id
                mlflow_run_name = run.data.tags.get('mlflow.runName')
                
                if swanlab.get_run() is None:
                    swanlab_run = swanlab.init(
                        project=self.project,
                        workspace=self.workspace,
                        experiment_name=mlflow_run_name,
                        mode=self.mode,
                        logdir=self.logdir,
                        description=run.data.tags.get('mlflow.note.content'),
                    )
                else:
                    swanlab_run = swanlab.get_run()
                
                swanlab_run.config.update(dict(
                    mlflow_run_id=run.info.run_id,
                    mlflow_run_name=mlflow_run_name,
                    mlflow_run_description=run.data.tags.get('mlflow.note.content'),
                    mlflow_run_params=run.data.params,
                    mlflow_run_tags={k: v for k, v in run.data.tags.items() if not k.startswith('mlflow')},
                ))

                # 预先获取所有 metric 历史数据，避免重复调用
                all_metric_histories = {}
                for key in run.data.metrics.keys():
                    history = client.get_metric_history(run_id, key)
                    all_metric_histories[key] = history

                # 计算总 metric 数量用于进度条
                total_metrics = sum(len(history) for history in all_metric_histories.values())

                progress = None
                task_id = None
                if Progress is not None:
                    progress = Progress(
                        TextColumn("[bold blue]{task.description}"),
                        BarColumn(bar_width=40),
                        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                        TimeElapsedColumn(),
                        TimeRemainingColumn(),
                    )
                    progress.start()
                    task_id = progress.add_task(mlflow_run_name or run_id[:8], total=total_metrics)

                index = 0
                for key, history in all_metric_histories.items():
                    for m in history:
                        swanlab_run.log({m.key: m.value}, step=m.step)
                        index += 1
                        if progress is not None and task_id is not None:
                            progress.update(task_id, advance=1)
                        # TODO: 等未来上传方案优化后解除延时
                        if index % 5 == 0:
                            time.sleep(1)

                if progress is not None:
                    progress.stop()
                
                swanlab_run.finish()
                
    def run(self, tracking_uri: str, experiment: int = None):
        swl.info("Start converting MLFlow Runs to SwanLab...")
        self.parse_mlflow_logs(
            tracking_uri=tracking_uri,
            experiment=experiment,
        )
        swl.info("Finished converting MLFlow Runs to SwanLab.")