#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/5/24 20:25
@File: pytest_env.py
@IDE: pycharm
@Description:
    测试swanlab/env.py中的环境变量和工具函数
"""
from swanlab.env import (
    SwanLabMode,
    get_mode,
    reset_env,
    MODE
)
import swanlab.env as E
from nanoid import generate
from tutils import SWANLAB_LOG_DIR, PACKAGE_PATH, TEMP_PATH
import pytest
import os


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    if MODE in os.environ:
        del os.environ[MODE]
    if E.ROOT in os.environ:
        del os.environ[E.ROOT]
    if E.PACKAGE in os.environ:
        del os.environ[E.PACKAGE]
    if E.DEV in os.environ:
        del os.environ[E.DEV]
    yield
    # 恢复原状
    if MODE in os.environ:
        del os.environ[MODE]


def use_strict_mode():
    """
    使用严格模式
    """
    get_mode({MODE: SwanLabMode.CLOUD.value})


def use_not_strict_mode():
    """
    使用非严格模式
    """
    get_mode({MODE: SwanLabMode.DISABLED.value})


class TestStrictMode:
    """
    同时测试get_mode和is_strict_mode函数
    """

    def test_strict_mode(self):
        """
        测试默认严格模式
        """
        assert E.is_strict_mode() is True

    def test_strict_mode_cloud(self):
        """
        测试严格模式
        """
        use_strict_mode()
        assert E.is_strict_mode() is True

    def test_strict_mode_disabled(self):
        """
        测试非严格模式
        """
        use_not_strict_mode()
        assert E.is_strict_mode() is False

    def test_strict_mode_local(self):
        """
        测试本地模式
        """
        use_strict_mode()
        assert E.is_strict_mode() is True

    def test_strict_mode_error(self):
        """
        测试重定向模式
        """
        with pytest.raises(ValueError):
            get_mode({MODE: generate()})


class TestAssertExist:
    """
    测试assert_exist函数的正确性，由于assert_exist函数有5个参数，除了第一个参数以外其他参数会影响返回值，所以需要分情况测试
    非严格模式与严格模式隔离，因此可以分开测试
    """

    def test_exist_non_strict(self):
        """
        非严格模式下的文件存在断言测试
        """
        use_not_strict_mode()
        assert E.assert_exist(__file__) is True
        assert E.assert_exist(generate()) is False

    # ---------------------------------- 严格模式测试 ----------------------------------

    def test_exist_ok(self):
        """
        测试严格模式下，文件路径存在
        """
        assert E.assert_exist(__file__) is True
        assert E.assert_exist(__file__, target_type="file") is True
        with pytest.raises(NotADirectoryError):
            E.assert_exist(__file__, target_type="folder")
        with pytest.raises(IsADirectoryError):
            E.assert_exist(os.path.dirname(__file__), target_type="file")
        with pytest.raises(FileNotFoundError):
            E.assert_exist(generate())

    def test_exist_ra_false(self):
        """
        测试严格模式下，文件路径不存在，但是不抛出异常
        """
        assert E.assert_exist(generate(), ra=False) is False
        assert E.assert_exist(generate(), target_type="file", ra=False) is False
        assert E.assert_exist(generate(), target_type="folder", ra=False) is False
        # 文件路径存在，但是设置了ra=False和target_type
        with pytest.raises(NotADirectoryError):
            E.assert_exist(__file__, target_type="folder", ra=False)
        with pytest.raises(IsADirectoryError):
            E.assert_exist(os.path.dirname(__file__), target_type="file", ra=False)


def test_db_path():
    """
    测试数据库路径
    """
    os.environ[E.ROOT] = SWANLAB_LOG_DIR
    assert E.get_db_path() == os.path.join(SWANLAB_LOG_DIR, "runs.swanlab")


def test_env_reset():
    os.environ[E.ROOT] = SWANLAB_LOG_DIR
    E.init_env()
    assert len(E._env) > 0
    reset_env()
    assert len(E._env) == 0


def test_is_dev():
    # 默认为False
    assert E.is_dev() is False
    reset_env()
    # 设置为True
    os.environ[E.DEV] = "TRUE"
    assert E.is_dev() is True
    reset_env()
    # 大小写敏感
    os.environ[E.DEV] = "true"
    assert E.is_dev() is False


class TestGetServer:
    """
    测试由环境变量设置服务器地址
    """

    def test_default_host(self):
        """
        测试默认地址
        """
        assert E.get_server_host() == "127.0.0.1"

    def test_use_env_host(self):
        """
        测试使用环境变量
        """
        os.environ[E.HOST] = "127.0.0.2"
        assert E.get_server_host() == "127.0.0.2"
        reset_env()
        os.environ[E.HOST] = generate()
        with pytest.raises(ValueError):
            E.get_server_host()

    def test_default_port(self):
        """
        测试默认端口
        """
        assert E.get_server_port() == 5092

    def test_use_env_port(self):
        """
        测试使用环境变量
        """
        os.environ[E.PORT] = "5093"
        assert E.get_server_port() == 5093
        reset_env()
        os.environ[E.PORT] = generate()
        with pytest.raises(ValueError):
            E.get_server_port()


class TestSwanLogDir:

    def test_default(self):
        """
        默认情况下日志目录应该存放在当前运行命令时的swanlog目录
        """
        use_not_strict_mode()
        assert E.get_swanlog_dir() == os.path.join(os.getcwd(), "swanlog")

    def test_use_env(self):
        """
        测试使用环境变量
        """
        os.environ[E.ROOT] = SWANLAB_LOG_DIR  # 在启动时已保证此目录存在
        assert E.get_swanlog_dir() == SWANLAB_LOG_DIR
        # 必须是一个绝对路径
        reset_env()
        del os.environ[E.ROOT]
        os.environ[E.ROOT] = generate()
        with pytest.raises(ValueError):
            E.get_swanlog_dir()

    def test_notfound_strict(self):
        """
        测试严格模式下目录不存在
        """
        use_strict_mode()
        with pytest.raises(FileNotFoundError):
            E.get_swanlog_dir()


class TestPackagePath:
    def test_default(self):
        """
        测试默认的package路径
        """
        project_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        package_path = os.path.join(project_dir, "swanlab", "package.json")
        assert E.get_package_path() == package_path

    def test_use_env(self):
        """
        非开发模式下使用环境变量
        """
        project_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        package_path = os.path.join(project_dir, "swanlab", "package.json")
        os.environ[E.PACKAGE] = generate()
        assert E.get_package_path() != os.environ[E.PACKAGE]
        assert E.get_package_path() == package_path

    def test_use_env_dev(self):
        """
        开发模式下使用环境变量
        """
        os.environ[E.DEV] = "TRUE"
        os.environ[E.PACKAGE] = generate()
        assert E.get_package_path() == os.environ[E.PACKAGE]
        # 需要注意的是测试时已经引入了swanlab包，而其默认执行了get_package_path函数
        from swanlab.package import package_path
        assert E.get_package_path() != package_path
        assert package_path == PACKAGE_PATH


class TestGetSwanLabFolder:
    """
    测试获取.swanlab目录路径
    """

    def test_use_env(self):
        os.environ[E.HOME] = TEMP_PATH
        assert TEMP_PATH not in E.get_swanlab_folder()

    def test_use_env_dev(self):
        os.environ[E.DEV] = "TRUE"
        os.environ[E.HOME] = TEMP_PATH
        assert E.get_swanlab_folder() == os.path.join(TEMP_PATH, ".swanlab")
