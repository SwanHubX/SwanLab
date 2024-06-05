#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_sdk.py
@IDE: pycharm
@Description:
    测试sdk的一些api
"""
from swanlab.env import (
    get_swanlab_folder,
    MODE,
    ROOT,
    HOME
)
import tutils as T
import swanlab.data.sdk as S
import swanlab.error as E
from swanlab.log import swanlog
from swanlab.data.run import get_run
from nanoid import generate
import pytest
import os


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    swanlog.disable_log()
    yield
    # 恢复原状
    swanlog.enable_log()
    if get_run() is not None:
        get_run().finish()
    os.environ[ROOT] = T.SWANLAB_LOG_DIR
    if HOME in os.environ:
        del os.environ[HOME]


class TestInitMode:
    """
    测试init时函数的mode参数设置行为
    """

    def test_init_disabled(self):
        run = S.init(mode="disabled", logdir=generate())
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})  # 不会报错
        a = run.settings.run_dir
        assert not os.path.exists(a)
        assert get_run() is not None

    def test_init_local(self):
        run = S.init(mode="local")
        assert os.environ[MODE] == "local"
        run.log({"TestInitMode": 1})  # 不会报错
        assert get_run() is not None

    def test_init_cloud(self):
        run = S.init(mode="cloud")
        assert os.environ[MODE] == "cloud"
        run.log({"TestInitMode": 1})  # 不会报错
        assert get_run() is not None

    def test_init_error(self):
        with pytest.raises(ValueError):
            S.init(mode="123456")  # noqa
        assert get_run() is None

    # ---------------------------------- 测试环境变量输入 ----------------------------------

    def test_init_disabled_env(self):
        os.environ[MODE] = "disabled"
        run = S.init()
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})
        a = run.settings.run_dir
        assert not os.path.exists(a)
        assert get_run() is not None

    def test_init_local_env(self):
        os.environ[MODE] = "local"
        run = S.init()
        assert os.environ[MODE] == "local"
        run.log({"TestInitMode": 1})

    def test_init_cloud_env(self):
        os.environ[MODE] = "cloud"
        run = S.init()
        assert os.environ[MODE] == "cloud"
        run.log({"TestInitMode": 1})

    # ---------------------------------- 环境变量和mode的情况 ----------------------------------

    def test_init_disabled_env_mode(self):
        os.environ[MODE] = "local"
        run = S.init(mode="disabled")
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})
        a = run.settings.run_dir
        assert not os.path.exists(a)
        assert get_run() is not None


class TestInitProject:
    """
    测试init时函数的project参数设置行为
    """

    def test_init_project_none(self):
        """
        设置project为None
        """
        run = S.init(project=None, mode="disabled")
        assert run.project_name == os.path.basename(os.getcwd())

    def test_init_project(self):
        """
        设置project为字符串
        """
        project = "test_project"
        run = S.init(project=project, mode="disabled")
        assert run.project_name == project


class TestInitLogdir:
    """
    测试init时函数的logdir参数设置行为
    """

    def test_init_logdir_disabled(self):
        """
        disabled模式下设置logdir不生效，采用的是环境变量的设置
        """
        logdir = generate()
        run = S.init(logdir=logdir, mode="disabled")
        assert run.settings.swanlog_dir != logdir
        assert run.settings.swanlog_dir == os.environ[ROOT]
        run.finish()
        del os.environ[ROOT]
        run = S.init(logdir=logdir, mode="disabled")
        assert run.settings.swanlog_dir != logdir
        assert run.settings.swanlog_dir == os.path.join(os.getcwd(), "swanlog")

    def test_init_logdir_enabled(self):
        """
        其他模式下设置logdir生效
        """
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        run = S.init(logdir=logdir, mode="local")
        assert run.settings.swanlog_dir == logdir
        run.finish()
        del os.environ[ROOT]
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        run = S.init(logdir=logdir, mode="local")
        assert run.settings.swanlog_dir == logdir

    def test_init_logdir_env(self):
        """
        通过环境变量设置logdir
        """
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[ROOT] = logdir
        run = S.init(mode="local")
        assert run.settings.swanlog_dir == logdir
        run.finish()
        del os.environ[ROOT]
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[ROOT] = logdir
        run = S.init(mode="local")
        assert run.settings.swanlog_dir == logdir


class TestLogin:
    """
    测试通过sdk封装的login函数登录
    不填apikey的部分不太好测
    """

    def test_use_home_key(self, monkeypatch):
        """
        使用家目录下的key，不需要输入
        如果家目录下的key获取失败，会使用getpass.getpass要求用户输入，作为测试，使用monkeypatch替换getpass.getpass
        """
        os.environ[HOME] = T.TEMP_PATH
        monkeypatch.setattr("getpass.getpass", T.get_password)
        S.login()
        # 默认保存Key
        assert os.path.exists(os.path.join(get_swanlab_folder(), ".netrc"))

    def test_use_input_key(self, monkeypatch):
        """
        使用输入的key
        """
        key = generate()
        with pytest.raises(E.ValidationError):
            S.login(api_key=key)
        key = T.KEY
        S.login(api_key=key)
