"""
@author: cunyue
@file: sync.py
@time: 2025/6/8 15:32
@description: 测试同步（仅测试本地解析）
"""

import os.path
import random
from typing import List

import numpy as np
import pytest
from nanoid import generate

import swanlab
from swanlab import sync
from swanlab.data.porter import DataPorter
from swanlab.toolkit import MetricInfo


def test_sync():
    """
    模拟一个本地实验，保存在本地，然后解析它
    """
    # ---------------------- 预备一些容器，存储每个产出，为后续验证做准备 ----------------
    record_metrics: List[MetricInfo] = []
    record_scalars: List[float] = []
    record_images_count: int = 0
    record_logs: List[str] = []
    # ---------------------------------- 开启实验 ----------------------------------
    config = {
        "epochs": 10000,
        "learning_rate": 0.01,
        "offset": random.random() / 5,
    }
    description = generate()
    tags = [generate(), generate()]
    # project_name = generate()
    project_name = "12345"
    # experiment_name = generate()
    experiment_name = "67890"
    run = swanlab.init(
        project=project_name,
        experiment_name=experiment_name,
        mode="offline",
        config=config,
        description=description,
        tags=tags,
    )
    for epoch in range(1, swanlab.config.epochs):
        acc = 1 - 2**-epoch - random.random() / epoch - swanlab.config.offset
        loss = 2**-epoch + random.random() / epoch + swanlab.config.offset
        ll = swanlab.log(
            {
                "accuracy": acc,
                "loss": loss,
                "image": swanlab.Image(np.random.random((3, 3, 3))),
            },
            step=epoch,
        )
        log = f"epoch={epoch}, accuracy={acc}, loss={loss}"
        print(log)
        record_scalars.extend([acc, loss])
        record_images_count += 1
        record_logs.append(log)
        record_metrics.extend([x for _, x in ll.items()])
    swanlab.finish()
    # ---------------------------------------------------------------------------
    run_dir = run.public.run_dir
    backup_file = run.public.backup_file
    # 文件夹存在
    assert os.path.exists(run_dir) is True
    # 解析日志文件成功（未登录）
    with pytest.raises(AssertionError) as e:
        sync(run_dir.__str__())
    assert e.value.args[0] == "Please log in first before using sync."
    assert os.path.isfile(backup_file), "Backup file does not exist after sync"
    # ---------------------- 验证所有的日志都存在 ----------------------------------
    # 1. 解析日志文件
    with DataPorter().open_for_sync(run_dir) as porter:
        porter.parse()
    header, project, experiment, logs, runtime, columns, scalars, medias, footer = (
        porter._header,
        porter._project,
        porter._experiment,
        porter._logs,
        porter._runtime,
        porter._columns,
        porter._scalars,
        porter._medias,
        porter._footer,
    )
    # 2.1 验证日志头部
    assert header.backup_type == "DEFAULT"
    # 2.2 验证项目内容
    assert project.name == project_name
    assert project.workspace is None
    assert project.public is None
    # 2.3 验证实验内容
    # 实验名称不为空，因为 swanlab 自动设置
    assert experiment.name is not None
    # 默认为空字符串
    assert experiment.description == description
    # 实验 tags
    assert experiment.tags == tags
    # 2.4 验证日志输出
    backup_logs = [log.message for log in logs]
    for record_log in record_logs:
        # 每一个 record_logs 都能在 logs 中找到
        # 注意判断顺序不能颠倒，因为 swanlab 本身会输出一些日志，这也会被记录到 logs 中
        assert record_log in backup_logs, "Log not found in backup logs: " + record_log
    # 2.5 验证运行时输出
    runtime_info = runtime.to_file_model(os.path.join(run_dir.__str__(), "files"))
    assert runtime_info.conda is None, "Not using conda, should be None"  # 此测试中没有开启 conda 检查
    assert isinstance(runtime_info.requirements, str), "Requirements should be a string"  # 运行时依赖检测成功
    assert isinstance(runtime_info.metadata, dict), "Metadata should be a dictionary"  # 系统元信息检测成功
    assert isinstance(runtime_info.config, dict), "Config should be a dictionary"  # 系统配置检测成功
    for key in runtime_info.config:
        # 验证配置项
        assert key in config, f"Config key {key} not found in original config"
        assert runtime_info.config[key]['value'] == config[key], f"Config value for {key} does not match original value"
    # 2.6 验证指标是否都存在
    assert len(scalars) + len(medias) == len(record_metrics), "Total metrics count does not match"
    # 验证标量类型
    backup_scalars = [
        metric.metric['data'] for metric in record_metrics if metric.column_info.chart_type.value.column_type == 'FLOAT'
    ]
    assert len(backup_scalars) == len(scalars), "Total scalars count does not match"
    for scalar in backup_scalars:
        # 验证标量值
        assert scalar in record_scalars, f"Scalar {scalar} not found in original scalars"
    # 验证媒体类型，这里简单一点，因为很难判断比如媒体类型的内容是否正确，所以只验证指标的数量和类型
    backup_images = [metric for metric in record_metrics if metric.column_info.chart_type.value.column_type == 'IMAGE']
    assert len(backup_images) == record_images_count, "Total images count does not match"
