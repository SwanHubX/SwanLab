"""
@author: cunyue
@file: test_client.py
@time: 2025/6/16 15:51
@description: 测试客户端功能
"""

import os

import nanoid
import pytest
import requests_mock

from swanlab.core_python import create_client, Client, CosClient, reset_client
from swanlab.core_python.auth import login_by_key
from swanlab.package import get_host_api
from swanlab.toolkit import MediaBuffer
from tutils import is_skip_cloud_test, TEMP_PATH, API_KEY
from tutils.setup import *


# ---------------------------------- mock 请求工具函数 ----------------------------------


class UseMocker(requests_mock.Mocker):
    """
    使用request_mock库进行mock测试，由于现在绝大部分请求都在get_host_api上，所以封装一层
    """

    def __init__(self, base_url: str = None):
        super().__init__()
        base_url = base_url or get_host_api()
        self.base_url = base_url

    def get(self, router, *args, **kwargs):
        return super().get(*(self.base_url + router, *args), **kwargs)

    def post(self, router, *args, **kwargs):
        return super().post(*(self.base_url + router, *args), **kwargs)

    def put(self, router, *args, **kwargs):
        return super().put(*(self.base_url + router, *args), **kwargs)

    def patch(self, router, *args, **kwargs):
        return super().patch(*(self.base_url + router, *args), **kwargs)

    def delete(self, router, *args, **kwargs):
        return super().delete(*(self.base_url + router, *args), **kwargs)


def test_use_mocker():
    with UseMocker() as m:
        m.post("/tmp", text="mock")
        import requests
        from swanlab.package import get_host_api

        resp = requests.post(get_host_api() + "/tmp")
        assert resp.text == "mock"


# ---------------------------------- 测试客户端 ----------------------------------


def test_decode_response():
    with UseMocker() as mocker:
        mocker.post("/json", json={"test": "test"})
        mocker.post("/text", text="test")
        with UseMockRunState() as run_state:
            client = run_state.client
            data, _ = client.post("/json")
            assert data == {"test": "test"}
            data, _ = client.post("/text")
            assert data == "test"


@pytest.mark.skipif(is_skip_cloud_test, reason="skip cloud test")
class TestCosSuite:
    http: Client = None
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    project_name = nanoid.generate(alphabet)
    experiment_name = nanoid.generate(alphabet)
    file_path = os.path.join(TEMP_PATH, nanoid.generate(alphabet))
    now_refresh_time = 1
    pre_refresh_time = CosClient.REFRESH_TIME

    @classmethod
    def setup_class(cls):
        CosClient.REFRESH_TIME = cls.now_refresh_time
        # 这里不测试保存token的功能
        login_info = login_by_key(API_KEY, save=False)
        cls.http = create_client(login_info)
        cls.http.mount_project(cls.project_name)
        cls.http.mount_exp(cls.experiment_name, ('#ffffff', '#ffffff'))
        # temp路径写一个文件上传
        with open(cls.file_path, "w") as f:
            f.write("test")

    @classmethod
    def teardown_class(cls):
        CosClient.REFRESH_TIME = cls.pre_refresh_time

    def test_cos_ok(self):
        assert self.http is not None
        assert self.http.cos is not None

    def test_cos_upload(self):
        # 新建一个文件对象
        buffer = MediaBuffer()
        buffer.write(b"test")
        buffer.file_name = "test"
        self.http.upload(buffer)
        # 为了开发方便，测试刷新功能关闭

        # # 开发版本设置的过期时间为3s，等待过期
        # time.sleep(3)
        # # 重新上传，测试刷新
        # assert self.http.cos.should_refresh is True
        # self.http.upload(buffer)


@pytest.mark.skipif(is_skip_cloud_test, reason="skip cloud test")
class TestExpSuite:
    """
    测试实验相关的功能
    """

    client: Client = None

    def setup_method(self):
        """
        初始化客户端对象
        """
        login_info = login_by_key(API_KEY, save=False)
        self.client = create_client(login_info)

    @staticmethod
    def teardown_method():
        """
        清理实验
        """
        reset_client()

    def test_create_exp_before_create_project(self):
        """
        实验创建前未创建项目，应该抛出异常
        """
        with pytest.raises(NotImplementedError) as e:
            self.client.mount_exp(nanoid.generate(), ('#ffffff', '#ffffff'))
        assert "Project not mounted, please call mount_project() first" == str(
            e.value
        ), "Expected ValueError when creating experiment without project"

    def test_create_exp_with_cuid(self):
        """
        测试创建实验时使用cuid
        cuid 必须为21位字符串且为小写字母和数字组成
        """
        cuid = nanoid.generate(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", size=21)
        exp_name = nanoid.generate()
        self.client.mount_project(nanoid.generate())
        self.client.mount_exp(exp_name, ('#ffffff', '#ffffff'), cuid=cuid)
        assert self.client.exp_id == cuid, "Experiment ID should match the provided cuid"

    def test_exp_must_exist(self):
        """
        测试实验在 cuid 传递时实验必须存在
        """
        proj_name = nanoid.generate()
        exp_name = nanoid.generate()
        self.client.mount_project(proj_name)
        cuid = nanoid.generate(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", size=21)
        with pytest.raises(AssertionError):
            # 如果传递 must_exist=True，则 cuid 必须传递
            self.client.mount_exp(exp_name, ('#ffffff', '#ffffff'), must_exist=True)
        with pytest.raises(RuntimeError):
            self.client.mount_exp(exp_name, ('#ffffff', '#ffffff'), cuid=cuid, must_exist=True)

    def test_exp_is_exist(self):
        """
        测试实验是否存在
        """
        proj_name = nanoid.generate()
        exp_name = "test"
        cuid = nanoid.generate(alphabet="abcdefghijklmnopqrstuvwxyz0123456789", size=21)
        self.client.mount_project(proj_name)
        new = self.client.mount_exp(exp_name, ('#ffffff', '#ffffff'), cuid=cuid)
        assert new is True, "Experiment should not exist before creation"
        new = self.client.mount_exp(exp_name, ('#ffffff', '#ffffff'), cuid=cuid)
        assert new is False, "Experiment should exist after creation"
