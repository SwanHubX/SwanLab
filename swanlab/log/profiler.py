"""
@author: Puiching-Memory
@file: profiler.py
@time: 2025/6/8 17:47
@description: 保存模型 profiler 日志
"""


def trace_handler(save_dir: str):
    """
    trace_handler 是一个回调函数，用于处理 torch.profiler 的 trace 信息，并将其保存到文件中

    examples
    -------
    >>> activities = [torch.profiler.ProfilerActivity.CPU, torch.profiler.ProfilerActivity.CUDA]
    >>> with torch.profiler.profile(activities=activities,on_trace_ready=trace_handler()) as p:
    """
    from . import swanlog
    import os

    assert os.path.isdir(save_dir), RuntimeError(
        "Run directory not found. Please ensure the run directory is properly set."
    )

    def handler_fn(prof) -> None:
        saved_path = os.path.join(save_dir, 'trace.json')
        if os.path.exists(saved_path):
            swanlog.warning(f"{saved_path} already exists, will be overwritten")
            os.remove(f"{saved_path}")
        else:
            swanlog.info(f"torch.profiler trace is saved to {saved_path}")

        prof.export_chrome_trace(f"{saved_path}")

    return handler_fn
