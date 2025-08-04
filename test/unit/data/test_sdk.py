#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: test_sdk.py
@IDE: pycharm
@Description:
    测试sdk的一些api
"""
import os
import random

import platformdirs
import pytest
from nanoid import generate

import swanlab.data.sdk as S
import swanlab.data.utils
import swanlab.error as Err
import tutils as T
from swanlab import Settings
from swanlab.core_python import reset_client, get_client
from swanlab.data.run import get_run
from swanlab.data.store import get_run_store
from swanlab.env import SwanLabEnv, get_save_dir
from swanlab.log import swanlog


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    swanlog.disable_log()
    run = get_run()
    if run is not None:
        run.finish()
    yield
    run = get_run()
    if run is not None:
        run.finish()
    swanlog.enable_log()


MODE = SwanLabEnv.MODE.value


class TestInitModeFunc:

    def test_init_error_mode(self):
        """
        初始化时mode参数错误
        """
        with pytest.raises(ValueError):
            swanlab.data.utils._init_mode("123456")  # noqa

    @pytest.mark.parametrize("mode", ["disabled", "local", "cloud"])
    def test_init_mode(self, mode, monkeypatch):
        """
        初始化时mode参数正确
        """
        if mode == 'cloud':
            mode = 'local'
            monkeypatch.setattr("builtins.input", lambda _: "3")
        swanlab.data.utils._init_mode(mode)
        assert os.environ[MODE] == mode
        del os.environ[MODE]
        # # 大写
        # S._init_mode(mode.upper())
        # assert os.environ[MODE] == mode

    @pytest.mark.parametrize("mode", ["disabled", "local", "cloud"])
    def test_overwrite_mode(self, mode, monkeypatch):
        """
        初始化时mode参数正确，覆盖环境变量
        """
        if mode == 'cloud':
            mode = 'local'
            monkeypatch.setattr("builtins.input", lambda _: "3")
        os.environ[MODE] = "123456"
        swanlab.data.utils._init_mode(mode)
        assert os.environ[MODE] == mode

    def test_no_api_key_to_cloud(self, monkeypatch):
        """
        初始化时mode为cloud，但是没有设置apikey
        """
        if SwanLabEnv.API_KEY.value in os.environ:
            del os.environ[SwanLabEnv.API_KEY.value]
        monkeypatch.setattr("builtins.input", lambda _: "3")
        mode, login_info = swanlab.data.utils._init_mode("cloud")
        assert mode == "offline"
        assert login_info is None

    @pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
    def test_init_cloud_with_no_api_key(self, monkeypatch):
        """
        初始化时mode为cloud，但是没有设置apikey
        """
        api_key = os.environ[SwanLabEnv.API_KEY.value]
        del os.environ[SwanLabEnv.API_KEY.value]
        # 在测试时默认会在交互模式下
        # 接下来需要模拟终端连接，使用monkeypatch
        # 三种选择方式：
        # 1. 输入api key
        # 2. 创建账号
        # 3. 使用本地版

        # 选择第三种
        monkeypatch.setattr("builtins.input", lambda _: "3")
        mode, login_info = swanlab.data.utils._init_mode("cloud")
        assert mode == "offline"
        assert login_info is None

        # 选择第二种
        monkeypatch.setattr("builtins.input", lambda _: "2")
        monkeypatch.setattr("getpass.getpass", lambda _: api_key)
        mode, login_info = swanlab.data.utils._init_mode("cloud")
        assert mode == "cloud"
        assert login_info is not None

        # 登录后会保存key，测试时需要删除
        os.remove(os.path.join(get_save_dir(), ".netrc"))

        # 选择第一种
        monkeypatch.setattr("builtins.input", lambda _: "1")
        monkeypatch.setattr("getpass.getpass", lambda _: api_key)
        mode, login_info = swanlab.data.utils._init_mode("cloud")
        assert mode == "cloud"
        assert login_info is not None

        # 登录后会保存key，测试时需要删除
        os.remove(os.path.join(get_save_dir(), ".netrc"))

        # 选择其他的
        def create_other_input():
            first = True

            def oi():
                nonlocal first
                if first:
                    first = False
                    return "123456"
                else:
                    return "3"

            return oi

        other_input = create_other_input()
        monkeypatch.setattr("builtins.input", lambda _: other_input())
        mode, login_info = swanlab.data.utils._init_mode("cloud")
        assert mode == "offline"
        assert login_info is None


class TestInitMode:
    """
    测试init时函数的mode参数设置行为
    """

    def test_init_local(self):
        run = S.init(mode="local")
        assert os.environ[MODE] == "local"
        run.log({"TestInitMode": 1})  # 不会报错
        assert run.mode == "local"
        assert run.public.cloud.available is False
        assert get_run() is not None
        assert run.public.cloud.project_name is None

    def test_init_disabled(self):
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        run = S.init(mode="disabled", logdir=logdir)
        assert not os.path.exists(logdir)
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})  # 不会报错
        a = run.public.run_dir
        assert os.path.exists(a)
        assert get_run() is not None

    @pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
    def test_init_cloud(self):
        S.login(T.is_skip_cloud_test)
        run = S.init(mode="cloud")
        assert os.environ[MODE] == "cloud"
        run.log({"TestInitMode": 1})  # 不会报错
        assert get_run() is not None
        for key in run.public.json()['cloud']:
            assert run.public.json()['cloud'][key] is not None

    def test_init_error(self):
        with pytest.raises(ValueError):
            S.init(mode="123456")  # noqa
        assert get_run() is None

    @pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
    def test_init_multiple(self):
        # 先初始化cloud
        self.test_init_cloud()
        get_run().finish()
        # 再初始化local
        self.test_init_local()
        get_run().finish()
        # 再初始化disabled
        self.test_init_disabled()

    # ---------------------------------- 测试环境变量输入 ----------------------------------

    def test_init_disabled_env(self):
        os.environ[MODE] = "disabled"
        run = S.init()
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})
        a = run.public.run_dir
        assert os.path.exists(a)
        assert get_run() is not None

    def test_init_local_env(self):
        os.environ[MODE] = "local"
        run = S.init()
        assert os.environ[MODE] == "local"
        run.log({"TestInitMode": 1})

    @pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
    def test_init_cloud_env(self):
        os.environ[MODE] = "cloud"
        S.login(T.is_skip_cloud_test)
        run = S.init()
        assert os.environ[MODE] == "cloud"
        run.log({"TestInitMode": 1})

    # ---------------------------------- 环境变量和mode的情况 ----------------------------------

    def test_init_disabled_env_mode(self):
        os.environ[MODE] = "local"
        run = S.init(mode="disabled")
        assert os.environ[MODE] == "disabled"
        run.log({"TestInitMode": 1})
        a = run.public.run_dir
        assert os.path.exists(a)
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
        assert run.public.project_name == os.path.basename(os.getcwd())

    def test_init_project(self):
        """
        设置project为字符串
        """
        project = "test_project"
        run = S.init(project=project, mode="disabled")
        assert run.public.project_name == project


def test_reinit():
    """
    测试重新初始化
    """
    run = S.init(mode="disabled")
    run_again = S.init(mode="local")
    assert run == run_again
    reinit = S.init(mode="disabled", reinit=True)
    assert run != reinit
    assert run.running is False
    assert run.success is True
    assert run.crashed is False


LOG_DIR = SwanLabEnv.SWANLOG_FOLDER.value


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
        assert run.public.swanlog_dir != logdir
        assert run.public.swanlog_dir == os.environ[LOG_DIR]
        run.finish()
        del os.environ[LOG_DIR]
        run = S.init(logdir=logdir, mode="disabled")
        assert run.public.swanlog_dir != logdir
        assert run.public.swanlog_dir == platformdirs.user_cache_dir(appname="swanlab", appauthor="SwanHubX")
        os.path.exists(run.public.swanlog_dir)

    def test_init_logdir_enabled(self):
        """
        其他模式下设置logdir生效
        """
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        run = S.init(logdir=logdir, mode="local")
        assert run.public.swanlog_dir == logdir
        run.finish()
        del os.environ[LOG_DIR]
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        run = S.init(logdir=logdir, mode="local")
        assert run.public.swanlog_dir == logdir

    def test_init_logdir_env(self):
        """
        通过环境变量设置logdir
        """
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[LOG_DIR] = logdir
        run = S.init(mode="local")
        assert run.public.swanlog_dir == logdir
        run.finish()
        del os.environ[LOG_DIR]
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[LOG_DIR] = logdir
        run = S.init(mode="local")
        assert run.public.swanlog_dir == logdir

    def test_init_logdir_offline(self):
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[LOG_DIR] = logdir
        run = S.init(mode="offline")
        assert run.public.swanlog_dir == logdir
        run.finish()
        del os.environ[LOG_DIR]
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[LOG_DIR] = logdir
        run = S.init(mode="offline")
        assert run.public.swanlog_dir == logdir

    @pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
    def test_init_logdir_settings_backup(self):
        logdir = os.path.join(T.TEMP_PATH, generate()).__str__()
        os.environ[LOG_DIR] = logdir
        # cloud模式下backup功能默认开启，此时依然生成文件夹
        run = S.init(mode="cloud")
        assert run.public.swanlog_dir == logdir
        os.path.exists(logdir)
        run.finish()
        # 如果通过 settings 设置为 False，则不会生成文件夹
        del os.environ[LOG_DIR]
        run = S.init(mode="cloud", settings=Settings(backup=False))
        assert run.public.swanlog_dir == platformdirs.user_cache_dir(appname="swanlab", appauthor="SwanHubX")
        assert run.public.swanlog_dir != logdir


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
class TestLogin:
    """
    测试通过sdk封装的login函数登录
    不填apikey的部分不太好测
    """

    @staticmethod
    def get_password(prompt: str):
        # 如果是第一次登录，使用错误的key，会提示重新输入
        if "Paste" in prompt:
            return generate()
        else:
            return T.is_skip_cloud_test

    def test_use_home_key(self, monkeypatch):
        """
        使用家目录下的key，不需要输入
        如果家目录下的key获取失败，会使用getpass.getpass要求用户输入，作为测试，使用monkeypatch替换getpass.getpass
        """
        os.environ[LOG_DIR] = T.TEMP_PATH
        monkeypatch.setattr("getpass.getpass", self.get_password)
        S.login()
        # 默认不保存Key
        assert not os.path.exists(os.path.join(get_save_dir(), ".netrc"))

    def test_use_home_key_save(self, monkeypatch):
        """
        使用家目录下的key，不需要输入
        如果家目录下的key获取失败，会使用getpass.getpass要求用户输入，作为测试，使用monkeypatch替换getpass.getpass
        """
        os.environ[LOG_DIR] = T.TEMP_PATH
        monkeypatch.setattr("getpass.getpass", self.get_password)
        S.login(save=True)
        # 保存Key
        assert os.path.exists(os.path.join(get_save_dir(), ".netrc"))

    def test_use_input_key(self, monkeypatch):
        """
        使用输入的key
        """
        key = generate(size=21, alphabet="12345")
        with pytest.raises(Err.ValidationError):
            S.login(api_key=key)
        key = T.API_KEY
        S.login(api_key=key)

    def test_use_env_key(self, monkeypatch):
        """
        测试code登录，使用环境变量key
        """

        def _():
            raise RuntimeError("this function should not be called")

        monkeypatch.setattr("getpass.getpass", _)
        os.environ[SwanLabEnv.API_KEY.value] = generate(size=21, alphabet="12345")
        with pytest.raises(Err.ValidationError):
            S.login()
        os.environ[SwanLabEnv.API_KEY.value] = T.API_KEY
        S.login()

    def test_pass_host(self):
        """
        传入host参数
        """
        del os.environ[SwanLabEnv.WEB_HOST.value]
        del os.environ[SwanLabEnv.API_HOST.value]
        os.environ[SwanLabEnv.API_KEY.value] = T.API_KEY
        S.login(host=T.API_HOST.rstrip("/api"))
        assert os.environ[SwanLabEnv.API_HOST.value] == T.API_HOST
        assert os.environ.get(SwanLabEnv.WEB_HOST.value) == T.API_HOST.rstrip("/api")
        S.login(host=T.API_HOST.rstrip("/api"), web_host=T.WEB_HOST)
        assert os.environ[SwanLabEnv.API_HOST.value] == T.API_HOST
        assert os.environ[SwanLabEnv.WEB_HOST.value] == T.WEB_HOST


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
class TestInitExpByEnv:
    """
    测试通过环境变量创建实验
    """

    def test_workspace(self):
        """
        测试通过环境变量创建实验
        """
        os.environ[SwanLabEnv.WORKSPACE.value] = generate()
        # 随便生成的一个workspace会报错
        with pytest.raises(ValueError) as e:
            S.init()
        assert str(e.value) == 'Space `{}` not found'.format(os.environ[SwanLabEnv.WORKSPACE.value])
        del os.environ[SwanLabEnv.WORKSPACE.value]
        reset_client()

    def test_exp_name(self):
        """
        测试通过环境变量创建实验
        """
        exp_name = generate()
        os.environ[SwanLabEnv.EXP_NAME.value] = exp_name
        run = S.init()
        assert run.public.cloud.experiment_name == exp_name
        del os.environ[SwanLabEnv.EXP_NAME.value]

    def test_proj_name(self):
        """
        测试通过环境变量创建实验
        """
        proj_name = generate()
        os.environ[SwanLabEnv.PROJ_NAME.value] = proj_name
        run = S.init()
        assert run.public.cloud.project_name == proj_name
        del os.environ[SwanLabEnv.PROJ_NAME.value]

    def test_env_params_priority(self):
        """
        环境变量的优先级低于函数参数
        """
        exp_name = generate()
        os.environ[SwanLabEnv.EXP_NAME.value] = exp_name
        param_exp_name = generate()
        run = S.init(experiment_name=param_exp_name)
        assert run.public.cloud.experiment_name == param_exp_name


@pytest.mark.timeout(10)
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_init_again_after_terminal_login(monkeypatch):
    """
    测试实验结束后再次开启新实验，此时不需要再次登录
    """
    api_key = os.environ[SwanLabEnv.API_KEY.value]
    del os.environ[SwanLabEnv.API_KEY.value]
    num = 0

    def choose_mode(_):
        nonlocal num
        if num == 0:
            num += 1
            return '1'
        else:
            raise RuntimeError("This function should not be called!")

    def input_api_key(_):
        nonlocal num
        if num == 1:
            num += 1
            return api_key
        else:
            raise RuntimeError("This function should not be called!")

    monkeypatch.setattr("builtins.input", choose_mode)
    monkeypatch.setattr("getpass.getpass", input_api_key)

    S.init()
    S.finish()

    def raise_error(_):
        raise RuntimeError("This function should not be called!")

    monkeypatch.setattr("builtins.input", raise_error)
    S.init()


@pytest.mark.timeout(10)
@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_init_again_after_code_login(monkeypatch):
    """
    测试在代码登录的情况下实验结束后再次开启新实验，此时不需要再次登录
    """

    def raise_error(_):
        raise RuntimeError("This function should not be called!")

    monkeypatch.setattr("builtins.input", raise_error)
    api_key = os.environ[SwanLabEnv.API_KEY.value]
    assert not os.path.exists(os.environ[SwanLabEnv.SWANLAB_FOLDER.value]), "SwanLab folder should be empty before test"

    del os.environ[SwanLabEnv.API_KEY.value]
    S.login(api_key=api_key)
    assert not os.path.exists(os.environ[SwanLabEnv.SWANLAB_FOLDER.value]), "SwanLab folder should be empty after login"

    S.init()
    S.finish()
    with pytest.raises(ValueError):
        get_client()
    assert os.getenv(SwanLabEnv.API_KEY.value) == api_key
    S.init()
    assert not os.path.exists(
        os.environ[SwanLabEnv.SWANLAB_FOLDER.value]
    ), "SwanLab folder should be empty after re-init"


@pytest.mark.parametrize("config", [1, [], (), True])
def test_init_error_config(config):
    """
    测试传入不同的config参数时，是否抛出错误对应错误
    """
    with pytest.raises(TypeError):
        swanlab.init(config=config, mode="disabled")


def test_init_error_config_name():
    """
    测试传入 string 类型参数
    此时会抛出 ValueError，因为不是指定的配置文件名后缀
    """
    with pytest.raises(ValueError):
        swanlab.init(mode="disabled", config="test_config")


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
class TestResume:
    """
    测试resume功能
    """

    @pytest.mark.parametrize("resume", ["must", "allow"])
    @pytest.mark.parametrize("mode", ["local", "disabled", "offline"])
    def test_allow_resume_with_invalid_mode(self, mode, resume):
        """
        除了 cloud 以外的模式不允许 resume='allow'
        """
        with pytest.raises(RuntimeError) as e:
            S.init(resume='allow', mode=mode)
        assert str(e.value) == "You can only use resume in cloud mode."

    @pytest.mark.parametrize("resume", ["must", "allow"])
    def test_sync_step(self, resume):
        """
        测试 resume 时同步 step 参数
        """
        run = S.init()
        run.log({"loss": 0.3, "accuracy": 0.7}, step=1)
        run.finish()
        import time

        time.sleep(5)
        new_run = S.init(resume=resume, id=run.id, reinit=True)
        assert new_run.id == run.id, "New run ID should be the same as the previous one"
        # 上传相同 step 的指标，此时报错
        ll = new_run.log({"loss": 0.3, "accuracy": 0.7}, step=1)
        assert ll["loss"].is_error
        assert ll["accuracy"].is_error
        # 上传大于 step=1 的指标，此时不会报错
        ll = new_run.log({"loss": 0.5, "accuracy": 0.8})
        assert not ll["loss"].is_error, "Expected loss metric to be logged successfully after reinit"
        assert not ll["accuracy"].is_error, "Expected accuracy metric to be logged successfully after reinit"
        assert ll["loss"].metric_step == 2, "Expected loss metric step to be 2 after reinit"
        assert ll["accuracy"].metric_step == 2, "Expected accuracy metric step to be 2 after reinit"

    # ---------------------------------- never 部分 ----------------------------------
    def test_never_resume(self):
        """
        测试never模式下的resume功能
        """
        run = S.init(resume='never')
        run_store = get_run_store()
        assert run.id is not None, "Run ID should be generated"
        assert run_store.run_id == run.id, "Run store ID should match the run ID"
        run.log({"loss": 0.1, "accuracy": 0.9}, step=1)
        run.finish()
        # 再次初始化，应该是一个新的实验
        new_run = S.init(resume='never')
        assert new_run.id != run.id, "New run ID should be different from the previous one"

    def test_never_with_id(self):
        """
        never 模式不允许传递 id 参数
        """
        with pytest.raises(RuntimeError) as e:
            S.init(resume='never', id='test_id')
        assert str(e.value) == "You can't pass id when resume=never or resume=False."

    # ---------------------------------- allow 部分 ----------------------------------

    def test_allow_resume(self):
        """
        测试allow模式下的resume功能
        """
        run = S.init(resume='allow')
        run_store = get_run_store()
        assert run.id is not None, "Run ID should be generated"
        assert run_store.run_id == run.id, "Run store ID should match the run ID"
        run.log({"loss": 0.1, "accuracy": 0.9}, step=1)
        run.finish()
        import time

        time.sleep(5)
        # 再次初始化，应该是一个新的实验
        new_run = S.init(resume='allow', id=run.id)
        assert new_run.id == run.id, "New run ID should be the same as the previous one"
        # 校验 run_store ID 是否在云端实验连接中
        assert run.id in new_run.public.cloud.experiment_url, "New run should be linked to the previous run"

    def test_allow_with_id(self):
        """
        allow 模式允许传递 id 参数，此时根据此 id 初始化实验
        """
        id = "".join(random.choices("abcdefghijklmnopqrstuvwxyz0123456789", k=21))
        run = S.init(resume='allow', id=id)
        run_store = get_run_store()
        assert run.id == id, "Run ID should match the provided ID"
        assert run_store.run_id == id, "Run store ID should match the provided ID"
        # 再次初始化实验
        new_run = S.init(resume='allow', id=id)
        assert new_run.id == id, "New run ID should be the same as the provided ID"

    # ---------------------------------- must 部分 ----------------------------------

    def test_must_resume(self):
        """
        must 模式下之前的实验必须存在，否则会报错
        """
        run = S.init()
        run.finish()
        # 再次初始化，应该是一个新的实验
        new_run = S.init(resume='must', id=run.id)
        assert new_run.id == run.id, "New run ID should be the same as the previous one"
        # 校验 run_store ID 是否在云端实验连接中
        assert run.id in new_run.public.cloud.experiment_url, "New run should be linked to the previous run"

    def test_must_no_id(self):
        """
        在 must 模式下必须传递 id 参数，否则会报错
        """
        with pytest.raises(ValueError) as e:
            S.init(resume='must')
        assert str(e.value) == "You must pass id when resume=must."


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
class TestResumeWithEnv:
    @staticmethod
    def teardown_method():
        """
        清理环境变量
        """
        if SwanLabEnv.RESUME.value in os.environ:
            del os.environ[SwanLabEnv.RESUME.value]

    def test_resume_env(self):
        """
        测试通过环境变量设置 resume 参数
        """
        os.environ[SwanLabEnv.RESUME.value] = "allow"
        S.init()
        run_store = get_run_store()
        assert run_store.resume == "allow"

    def test_resume_env_passed(self):
        """
        显式传递 resume 参数时，环境变量中的 resume 参数不生效
        """
        os.environ[SwanLabEnv.RESUME.value] = "allow"
        S.init(resume="never")
        run_store = get_run_store()
        assert run_store.resume == "never"


@pytest.mark.skipif(T.is_skip_cloud_test, reason="skip cloud test")
def test_online():
    """
    测试在线模式
    """
    run = S.init(mode="online")
    assert run.mode == "cloud"
    run.finish()
