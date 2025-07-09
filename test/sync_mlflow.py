import mlflow
import random
import swanlab

swanlab.sync_mlflow()

mlflow.set_experiment("mlflow_sync_test")

with mlflow.start_run(run_name="test_run"):
    mlflow.log_param("learning_rate", 0.01)
    mlflow.log_params({"batch_size": 32, "epochs": 10})
    
    for epoch in range(10):
        acc = 1 - 2 ** -epoch - random.random() / epoch
        loss = 2 ** -epoch + random.random() / epoch
        mlflow.log_metric("accuracy", acc, step=epoch)
        mlflow.log_metric("loss", loss, step=epoch)
        
        mlflow.log_metrics({
            "precision": acc * 0.9,
            "recall": acc * 0.8
        }, step=epoch)