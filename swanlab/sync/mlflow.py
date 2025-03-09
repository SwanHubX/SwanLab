"""
MLflow synchronization module, used to synchronize MLflow operations to SwanLab.

Test script:
```python
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
"""

import swanlab

def _extract_args(args, kwargs, param_names):
    """
    从args和kwargs中提取参数值的通用函数
    
    Args:
        args: 位置参数元组
        kwargs: 关键字参数字典
        param_names: 参数名称列表
    
    Returns:
        tuple: 按param_names顺序返回提取的参数值
    """
    values = []
    for i, name in enumerate(param_names):
        if len(args) > i:
            values.append(args[i])
        else:
            values.append(kwargs.get(name, None))
    return tuple(values)


def sync_mlflow(mode:str="cloud", mlflow_tracking:bool=True):
    """
    同步MLflow到SwanLab，支持记录参数、指标和模型
    
    Args:
        mode: "cloud", "local" 或 "disabled"。详见 https://docs.swanlab.cn/api/py-init.html
        mlflow_tracking: 如果此参数设置为False，则不会将数据上传到MLflow服务器
    
    使用示例:
    ```python
    import mlflow
    import random
    import swanlab
    
    swanlab.sync_mlflow()
    
    mlflow.set_experiment("my_experiment")
    
    with mlflow.start_run(run_name="my_run"):
        # 记录参数
        mlflow.log_param("learning_rate", 0.01)
        mlflow.log_params({"batch_size": 32, "epochs": 10})
        
        # 记录指标
        for epoch in range(10):
            acc = 1 - 2 ** -epoch - random.random() / epoch
            loss = 2 ** -epoch + random.random() / epoch
            mlflow.log_metric("accuracy", acc, step=epoch)
            mlflow.log_metric("loss", loss, step=epoch)
    ```
    """
    try:
        import mlflow
        from mlflow.tracking import MlflowClient
    except ImportError:
        raise ImportError("请先安装MLflow，命令: `pip install mlflow`")
    
    # 保存原始函数
    original_start_run = mlflow.start_run
    original_end_run = mlflow.end_run
    original_log_param = mlflow.log_param
    original_log_params = mlflow.log_params
    original_log_metric = mlflow.log_metric
    original_log_metrics = mlflow.log_metrics
    original_set_experiment = mlflow.set_experiment
    
    def patched_set_experiment(experiment_name, *args, **kwargs):
        """拦截设置实验名称的调用"""
        # 如果SwanLab尚未初始化，则使用MLflow的实验名称初始化
        if swanlab.data.get_run() is None:
            swanlab.init(project=experiment_name, mode=mode)
        return original_set_experiment(experiment_name, *args, **kwargs)
    
    def patched_start_run(*args, **kwargs):
        """拦截开始运行的调用"""
        run_id, experiment_id, run_name, nested, tags = _extract_args(
            args, kwargs, ['run_id', 'experiment_id', 'run_name', 'nested', 'tags']
        )
        
        # 如果SwanLab尚未初始化，则初始化它
        if swanlab.data.get_run() is None:
            swanlab.init(experiment_name=run_name, mode=mode)
        elif run_name:  # 如果已初始化但提供了新的run_name
            # 更新配置以包含run_name
            swanlab.config.update({"mlflow_run_name": run_name})
        
        # 如果提供了tags，将其添加到SwanLab配置
        if tags:
            swanlab.config.update({"mlflow_tags": tags})
        
        if not mlflow_tracking:
            # 如果不需要MLflow跟踪，则设置为不跟踪模式
            # 注意：MLflow没有直接的离线模式，这里可能需要根据MLflow的API调整
            pass
        
        return original_start_run(*args, **kwargs)
    
    def patched_end_run(*args, **kwargs):
        """拦截结束运行的调用"""
        status = _extract_args(args, kwargs, ['status'])[0]
        
        # 结束SwanLab运行
        swanlab.finish()
        
        return original_end_run(*args, **kwargs)
    
    def patched_log_param(key, value, *args, **kwargs):
        """拦截记录单个参数的调用"""
        # 将参数添加到SwanLab配置
        swanlab.config.update({key: value})
        
        return original_log_param(key, value, *args, **kwargs)
    
    def patched_log_params(params, *args, **kwargs):
        """拦截记录多个参数的调用"""
        # 将参数添加到SwanLab配置
        swanlab.config.update(params)
        
        return original_log_params(params, *args, **kwargs)
    
    def patched_log_metric(key, value, step=None, *args, **kwargs):
        """拦截记录单个指标的调用"""
        # 记录到SwanLab
        swanlab.log({key: value}, step=step)
        
        return original_log_metric(key, value, step=step, *args, **kwargs)
    
    def patched_log_metrics(metrics, step=None, *args, **kwargs):
        """拦截记录多个指标的调用"""
        # 记录到SwanLab
        swanlab.log(metrics, step=step)
        
        return original_log_metrics(metrics, step=step, *args, **kwargs)
    
    # 应用monkey patch
    mlflow.start_run = patched_start_run
    mlflow.end_run = patched_end_run
    mlflow.log_param = patched_log_param
    mlflow.log_params = patched_log_params
    mlflow.log_metric = patched_log_metric
    mlflow.log_metrics = patched_log_metrics
    mlflow.set_experiment = patched_set_experiment
