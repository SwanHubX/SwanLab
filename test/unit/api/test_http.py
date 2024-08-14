#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/31 12:47
@File: pytest_cos.py
@IDE: pycharm
@Description:
    测试http的api
    开发环境下存储凭证过期时间为3s
"""
import os
import nanoid
from swanlab.api.http import create_http, HTTP, CosClient
from swanlab.api.auth.login import login_by_key
from swanlab.data.modules import MediaBuffer
from tutils import API_KEY, TEMP_PATH, is_skip_cloud_test
from tutils.setup import UseMocker, UseSetupHttp
import pytest
import responses
from responses import registries

alphabet = "abcdefghijklmnopqrstuvwxyz"


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
    http: HTTP = None
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
        cls.http = create_http(login_info)
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
