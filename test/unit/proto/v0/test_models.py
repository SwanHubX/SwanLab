"""
@author: cunyue
@file: test_models.py
@time: 2025/6/8 14:16
@description: 测试备份模型
"""

import pytest

from swanlab.data.backup import BaseModel, Log, Runtime
from swanlab.data.backup.writer import write_runtime_info
from swanlab.log.type import LogData, LogContent
from swanlab.toolkit import create_time
from tutils import TEMP_PATH


def test_base():
    base_model = BaseModel()
    assert base_model.to_record() == '{"model_type": "BaseModel", "data": {}}\n'
    # BaseModel 本身无法使用 from_record 方法进行反序列化
    with pytest.raises(ValueError):
        assert base_model.from_record('{"model_type": "BaseModel", "data": {}}\n') is not None
    with pytest.raises(ValueError):
        base_model.from_record('{"model_type": "WrongType", "data": {}\n')


@pytest.mark.parametrize(
    "log_type",
    [
        ["stdout", "INFO"],
        ["stderr", "WARN"],
    ],
)
def test_log(log_type):
    ct = create_time()
    log_data = LogData(type=log_type[0], contents=[LogContent(message="Test log message", create_time=ct, epoch=1)])
    logs = Log.from_log_data(log_data=log_data)
    assert len(logs) == 1
    log_model = logs[0].to_log_model()
    assert log_model['level'] == log_type[1]
    assert len(log_model['contents']) == 1
    log_content = log_model['contents'][0]
    assert log_content['message'] == "Test log message"
    assert log_content['epoch'] == 1
    assert log_content['create_time'] == ct


@pytest.mark.parametrize("conda", [None, "test conda"])
@pytest.mark.parametrize("requirements", [None, "test requirements"])
@pytest.mark.parametrize("metadata", [None, {"key": "value"}])
@pytest.mark.parametrize("config", [None, {"setting": "value"}])
def test_runtime(conda, requirements, metadata, config):
    from swanlab.toolkit import RuntimeInfo

    # 1. runtime info -> runtime
    runtime_info = RuntimeInfo(requirements, metadata, config, conda)
    #  模拟文件写入
    write_runtime_info(TEMP_PATH, runtime_info)
    runtime = Runtime.from_runtime_info(runtime_info)
    if conda is not None:
        assert runtime.conda_filename == runtime_info.conda.name
    else:
        assert runtime.conda_filename is None
    if requirements is not None:
        assert runtime.requirements_filename == runtime_info.requirements.name
    else:
        assert runtime.requirements_filename is None
    if metadata is not None:
        assert runtime.metadata_filename == runtime_info.metadata.name
    else:
        assert runtime.metadata_filename is None
    if config is not None:
        assert runtime.config_filename == runtime_info.config.name
    else:
        assert runtime.config_filename is None
    # 2. runtime -> file model
    file_model = runtime.to_file_model(TEMP_PATH)
    if conda is not None:
        file_model.conda = runtime_info.conda.name
    else:
        assert file_model.conda is None
    if requirements is not None:
        file_model.requirements = runtime_info.requirements.name
    else:
        assert file_model.requirements is None
    if metadata is not None:
        file_model.metadata = runtime_info.metadata.name
    else:
        assert file_model.metadata is None
    if config is not None:
        file_model.config = runtime_info.config.name
    else:
        assert file_model.config is None


# 其他类似 但是感觉没必要写
