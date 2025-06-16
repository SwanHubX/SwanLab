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
import responses
from responses import registries
from swankit.core import MediaBuffer

from swanlab.api.auth.login import login_by_key
from swanlab.core_python import create_client, Client, CosClient
from swanlab.package import get_host_api
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
        with UseSetupHttp() as http:
            data = http.post("/json")
            assert data == {"test": "test"}
            data = http.post("/text")
            assert data == "test"


@responses.activate(registry=registries.OrderedRegistry)
def test_retry():
    """
    测试重试机制
    """
    from swanlab.package import get_host_api

    url = get_host_api() + "/retry"
    rsp1 = responses.get(url, body="Error", status=500)
    rsp2 = responses.get(url, body="Error", status=500)
    rsp3 = responses.get(url, body="Error", status=500)
    rsp4 = responses.get(url, body="OK", status=200)
    with UseSetupHttp() as http:
        data = http.get("/retry")
        assert data == "OK"
        assert rsp1.call_count == 1
        assert rsp2.call_count == 1
        assert rsp3.call_count == 1
        assert rsp4.call_count == 1


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
