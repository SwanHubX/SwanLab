"""
------example.py------
from swanlab.converter import MLFlowConverter

mlf_converter = MLFlowConverter()
mlf_converter.run(tracking_uri="MLFLOW_TRACKING_URL", experiment=0)
"""
import swanlab
from swanlab.log import swanlog as swl

class MLFLowConverter:
    def __init__(
        self,
        project: str = None,
        workspace: str = None,
        cloud: bool = True,
        logdir: str = None,
        
        **kwargs,
    ):
        
        self.project = project
        self.workspace = workspace
        self.cloud = cloud
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
                        mode="cloud" if self.cloud else "local",
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
                
                for key in run.data.metrics.keys():
                    for m in client.get_metric_history(run_id, key):
                        swanlab_run.log({m.key: m.value}, step=m.step)
                
                swanlab_run.finish()
                
    def run(self, tracking_uri: str, experiment: int = None):
        swl.info("Start converting MLFlow Runs to SwanLab...")
        self.parse_mlflow_logs(
            tracking_uri=tracking_uri,
            experiment=experiment,
        )
        swl.info("Finished converting MLFlow Runs to SwanLab.")