"""
@author: cunyue
@file: test_cloud.py
@time: 2025/6/26 22:55
@description: 测试云端回调功能
"""

import os

import pytest

import swanlab
import tutils as T
from swanlab.core_python import reset_client, get_client
from swanlab.core_python.auth.providers.api_key import login_by_key
from swanlab.data.callbacker import CloudPyCallback
from swanlab.data.store import get_run_store
from swanlab.env import SwanLabEnv
from tutils.setup import mock_login_info, UseMockRunState


def test_create_http(monkeypatch):
    """
    测试登录信息存在时的行为
    """
    login_info = mock_login_info()
    if os.getenv(SwanLabEnv.API_KEY.value):
        del os.environ[SwanLabEnv.API_KEY.value]
    with UseMockRunState():
        reset_client()
        with pytest.raises(ValueError):
            get_client()

        # 默认情况下会报错，因为没有登录信息
        def raise_error(_):
            raise RuntimeError("Mocked error for testing")

        monkeypatch.setattr("getpass.getpass", raise_error)
        with pytest.raises(RuntimeError) as e:
            CloudPyCallback._create_client()
        assert str(e.value) == "Mocked error for testing"

        # 写入登录信息后不报错
        CloudPyCallback.login_info = login_info
        callback = CloudPyCallback()
        client = callback._create_client()
        assert client == get_client()
        assert CloudPyCallback.login_info is None


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_recreate_http(monkeypatch):
    """
    测试登录信息不存在时的行为
    1. 如果 CloudPyCallback.login_info 存在，则使用该登录信息创建客户端
    2. 如果环境变量中存在API_KEY，则使用该API_KEY登录
    3. 否则让进入交互，让用户输入API_KEY登录
    """
    api_key = os.getenv(SwanLabEnv.API_KEY.value)
    del os.environ[SwanLabEnv.API_KEY.value]

    # 输入 api-key 登录
    def input_api_key(_):
        return api_key

    def raise_error(_):
        raise RuntimeError("Should not be called")

    # 1. 如果 CloudPyCallback.login_info 存在，则使用该登录信息创建客户端
    monkeypatch.setattr("getpass.getpass", raise_error)
    login_info = login_by_key(api_key, save=False)
    with pytest.raises(RuntimeError) as e:
        CloudPyCallback._create_client()
    assert str(e.value) == "Should not be called"
    CloudPyCallback.login_info = login_info
    callback = CloudPyCallback()
    client = callback._create_client()
    assert CloudPyCallback.login_info is None
    assert client == get_client()
    callback.porter.close()
    reset_client()
    with pytest.raises(ValueError):
        get_client()
    # 2. 如果环境变量中存在API_KEY，则使用该API_KEY登录
    callback = CloudPyCallback()
    assert callback.login_info is None
    os.environ[SwanLabEnv.API_KEY.value] = api_key
    client2 = callback._create_client()
    assert client2 == get_client()
    assert client != client2
    callback.porter.close()
    reset_client()

    # 3. 否则让进入交互，让用户输入API_KEY登录
    monkeypatch.setattr("getpass.getpass", input_api_key)
    callback = CloudPyCallback()
    assert callback.login_info is None
    client3 = callback._create_client()
    assert client3 == get_client()
    assert client2 != client3 != client


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
class TestCloudResume:
    """
    测试云端恢复功能
    """

    @staticmethod
    def mock_run_state(workspace: str, project="test-cloud-resume"):
        run_store = get_run_store()
        run_store.project = project
        run_store.workspace = workspace
        run_store.visibility = True
        return run_store

    def test_init_config(self):
        """
        测试配置初始化
        """
        config = {
            "name": "test_cloud_resume",
            "description": "This is a test experiment for cloud resume.",
            "epochs": 10,
            "learning_rate": 0.001,
            "batch_size": 32,
        }
        run = swanlab.init(config=config, project="test-cloud-resume")
        run.finish()
        for key in run.config:
            assert run.config[key] == config[key]
        swanlab.login()
        client = get_client()
        # resume 实验
        with UseMockRunState(run_id=run.id, client=client):
            run_store = self.mock_run_state(client.username, run.public.project_name)
            run_store.resume = "must"
            callback = CloudPyCallback()
            callback.on_init()
            assert run_store.config is not None, "Config should be loaded when resuming."
            assert run_store.config == config, "Config should be loaded when resuming."
