#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/18 15:28
@File: test_env.py
@IDE: pycharm
@Description:
    测试swanlab.env模块
"""
import datetime
import os
import shutil
from pathlib import Path

import pytest

import swanlab
import swanlab.env as E
from swanlab.env import (
    SwanLabEnv,
    is_interactive,
    get_save_dir,
    remove_host_suffix,
    create_swanlog_dir,
    create_time,
    SwanLabMode,
    utc_time,
)


def test_utc_time():
    t = utc_time()
    assert isinstance(t, datetime.datetime)
    assert t.tzinfo == datetime.timezone.utc
    first_call = utc_time()
    second_call = utc_time()
    assert second_call >= first_call


def test_create_time():
    t = create_time()
    assert t.endswith("+00:00")
    d = datetime.datetime.fromisoformat(t)
    assert d.tzinfo == datetime.timezone.utc


def test_list_mode():
    ms = SwanLabMode.list()
    assert len(ms) == 4
    assert "disabled" in ms
    assert "cloud" in ms
    assert "offline" in ms
    assert "local" in ms


class TestGetFolder:
    """
    测试获取全局保存文件夹
    """

    def test_default(self):
        """
        默认情况
        """
        if SwanLabEnv.SWANLAB_FOLDER.value in os.environ:
            del os.environ[SwanLabEnv.SWANLAB_FOLDER.value]
        # 获取当前用户的主目录
        home = os.path.expanduser("~")
        folder = os.path.join(home, ".swanlab")
        if os.path.exists(folder):
            shutil.rmtree(folder)
        assert E.get_save_dir() == folder
        assert os.path.exists(folder)

    def test_env_abs(self, tmp_path):
        """
        设置了一个绝对路径的环境变量
        """
        path = tmp_path / "env"
        os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = path.__str__()
        assert Path(E.get_save_dir()) == path
        assert os.path.exists(path)

    def test_env_rel(self, tmp_path):
        """
        设置了一个相对路径的环境变量
        """
        abs_path = tmp_path / "env_rel"
        rel_path = os.path.relpath(abs_path.__str__(), os.getcwd())
        os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = rel_path
        assert not os.path.exists(abs_path)
        assert E.get_save_dir() == abs_path.__str__()
        assert os.path.exists(abs_path)

    def test_env_parent_not_exist(self, tmp_path):
        """
        父目录不存在
        """
        path = tmp_path / "a" / "b"
        os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = path.__str__()
        assert not os.path.exists(path)
        with pytest.raises(FileNotFoundError):
            E.get_save_dir()

    def test_env_not_a_folder(self, tmp_path):
        """
        文件夹不存在，但是文件存在
        """
        path = tmp_path / "a_file"
        with open(path, "w") as f:
            f.write("test")
        os.environ[SwanLabEnv.SWANLAB_FOLDER.value] = path.__str__()
        assert os.path.exists(path)
        with pytest.raises(NotADirectoryError):
            E.get_save_dir()


class TestGetSwanlog:

    def test_default(self):
        """
        默认情况
        """
        if SwanLabEnv.SWANLOG_FOLDER.value in os.environ:
            del os.environ[SwanLabEnv.SWANLOG_FOLDER.value]
        pwd = os.getcwd()
        folder = os.path.join(pwd, "swanlog")
        assert E.get_swanlog_dir() == folder
        assert not os.path.exists(folder)

    def test_env_abs(self, tmp_path):
        """
        设置了一个绝对路径的环境变量
        """
        path = tmp_path / "env_swanlog"
        os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = path.__str__()
        assert E.get_swanlog_dir() == path.__str__()
        assert not os.path.exists(path)

    def test_env_rel(self, tmp_path):
        """
        设置了一个相对路径的环境变量
        """
        abs_path = tmp_path / "env_rel"
        rel_path = os.path.relpath(abs_path.__str__(), os.getcwd())
        os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = rel_path.__str__()
        assert E.get_swanlog_dir() == abs_path.__str__()
        assert not os.path.exists(abs_path)

    def test_env_parent_not_exist(self, tmp_path):
        """
        父目录不存在
        """
        path = tmp_path / "a" / "b"
        os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = path.__str__()
        assert not os.path.exists(path)
        with pytest.raises(FileNotFoundError):
            E.get_swanlog_dir()

    def test_env_not_a_folder(self, tmp_path):
        """
        文件夹不存在，但是文件存在
        """
        path = tmp_path / "a_file"
        with open(path, "w") as f:
            f.write("test")
        os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = path.__str__()
        assert os.path.exists(path)
        with pytest.raises(NotADirectoryError):
            E.get_swanlog_dir()


class TestGetMode:
    """
    测试获取解析模式
    """

    def test_default(self):
        """
        默认情况
        """
        assert E.get_mode() == E.SwanLabMode.CLOUD.value

    def test_env(self):
        """
        设置了环境变量
        """
        os.environ[SwanLabEnv.MODE.value] = E.SwanLabMode.LOCAL.value
        assert E.get_mode() == E.SwanLabMode.LOCAL.value

    def test_unknown(self):
        """
        未知的模式
        """
        os.environ[SwanLabEnv.MODE.value] = "unknown"
        with pytest.raises(ValueError):
            E.get_mode()


def test_set_default():
    """
    测试获取默认的环境变量
    """
    del os.environ[SwanLabEnv.WEB_HOST.value]
    del os.environ[SwanLabEnv.API_HOST.value]
    del os.environ[SwanLabEnv.RUNTIME.value]
    swanlab.env.SwanLabEnv.set_default()
    assert swanlab.package.get_host_web() == "https://swanlab.cn"
    assert swanlab.package.get_host_api() == "https://api.swanlab.cn/api"
    assert os.getenv(SwanLabEnv.RUNTIME.value) == "user"


def test_set_by_netrc():
    """
    测试通过netrc文件设置环境变量
    """
    del os.environ[SwanLabEnv.WEB_HOST.value]
    del os.environ[SwanLabEnv.API_HOST.value]
    del os.environ[SwanLabEnv.RUNTIME.value]
    netrc_path = os.path.join(get_save_dir(), ".netrc")
    with open(netrc_path, "w") as f:
        f.write("machine https://example.ai\nlogin test\npassword 123")
    swanlab.env.SwanLabEnv.set_default()
    assert swanlab.package.get_host_web() == "https://example.ai"
    assert swanlab.package.get_host_api() == "https://example.ai/api"
    assert os.getenv(SwanLabEnv.RUNTIME.value) == "user"


def test_check():
    """
    测试检查环境变量
    """
    os.environ[SwanLabEnv.MODE.value] = "124345"
    with pytest.raises(ValueError):
        SwanLabEnv.check()
    os.environ[SwanLabEnv.RUNTIME.value] = "124"
    with pytest.raises(ValueError):
        SwanLabEnv.check()


def test_is_interactive():
    # 测试时默认返回true
    assert is_interactive() == True


# ---------------------------------- 测试移除 host 后缀 ----------------------------------
def test_removes_suffix_when_host_ends_with_suffix():
    host = "example.ai/api"
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == "example.ai"


def test_returns_original_host_when_host_does_not_end_with_suffix():
    host = "example.com"
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == "example.com"


def test_handles_empty_host_and_returns_empty():
    host = ""
    suffix = "/api"
    assert remove_host_suffix(host, suffix) == ""


def test_handles_empty_suffix_and_returns_original_host():
    host = "example.com/api"
    suffix = ""
    assert remove_host_suffix(host, suffix) == "example.com/api"


def test_handles_suffix_longer_than_host_and_returns_original_host():
    host = "api"
    suffix = "example.com/api"
    assert remove_host_suffix(host, suffix) == "api"


def test_handles_host_with_multiple_suffixes():
    host = "example.com/api/v1"
    suffix = "/api/v1"
    assert remove_host_suffix(host, suffix) == "example.com"


def test_handles_host_with_blank_suffix():
    host = "example.com  "
    suffix = "api"
    assert remove_host_suffix(host, suffix) == "example.com"


# 测试创建新目录
def test_create_new_directory(tmp_path):
    """测试创建新目录的功能"""
    test_dir = tmp_path / "new_swanlog"
    result = create_swanlog_dir(test_dir)

    assert result == test_dir
    assert os.path.exists(test_dir)
    assert os.path.isdir(test_dir)
    assert os.access(test_dir, os.W_OK)


# 测试目录已存在的情况
def test_existing_directory(tmp_path):
    """测试目录已存在的情况"""
    test_dir = tmp_path / "existing_swanlog"
    test_dir.mkdir()  # 预先创建目录

    result = create_swanlog_dir(test_dir)

    assert result == test_dir
    assert os.path.exists(test_dir)


# 测试创建 .gitignore 文件
def test_gitignore_creation(tmp_path):
    """测试在空目录中创建 .gitignore 文件"""
    test_dir = tmp_path / "empty_swanlog"
    test_dir.mkdir()

    create_swanlog_dir(test_dir)

    gitignore = test_dir / ".gitignore"
    assert gitignore.exists()
    with open(gitignore, "r", encoding="utf-8") as f:
        assert f.read() == "*"


# 测试非空目录不创建 .gitignore
def test_non_empty_directory(tmp_path):
    """测试非空目录不创建 .gitignore"""
    test_dir = tmp_path / "non_empty_swanlog"
    test_dir.mkdir()
    (test_dir / "existing_file.txt").touch()  # 创建文件使目录非空

    create_swanlog_dir(test_dir)

    gitignore = test_dir / ".gitignore"
    assert not gitignore.exists()  # 不应创建 .gitignore


# 测试默认参数行为
def test_default_parameter(monkeypatch, tmp_path):
    """测试使用默认参数"""

    # 使用 monkeypatch 模拟 get_swanlog_dir 函数
    def mock_get_swanlog_dir():
        return tmp_path / "default_swanlog"

    monkeypatch.setattr("swanlab.env.get_swanlog_dir", mock_get_swanlog_dir)

    result = create_swanlog_dir()

    expected_dir = tmp_path / "default_swanlog"
    assert result == expected_dir
    assert os.path.exists(expected_dir)


# 测试无写入权限的情况
def test_no_write_permission(tmp_path):
    """测试无写入权限的情况"""
    test_dir = tmp_path / "no_permission"
    test_dir.mkdir()

    # 修改目录权限为只读
    os.chmod(test_dir, 0o444)  # 只读权限

    with pytest.raises(IOError) as exc_info:
        create_swanlog_dir(test_dir)

    assert "no write permission" in str(exc_info.value)

    # 恢复权限以便临时目录可以清理
    os.chmod(test_dir, 0o755)


# 测试创建目录失败的情况
def test_directory_creation_failure(monkeypatch, tmp_path):
    """测试目录创建失败的情况"""
    test_dir = tmp_path / "failed_creation"

    # 使用 monkeypatch 模拟 os.makedirs 抛出异常
    def mock_makedirs(*_, **__):
        raise OSError("Mocked creation error")

    monkeypatch.setattr("os.makedirs", mock_makedirs)

    with pytest.raises(IOError) as exc_info:
        create_swanlog_dir(test_dir)

    assert "Failed to create or access logdir" in str(exc_info.value)
    assert "Mocked creation error" in str(exc_info.value)


# 测试各种路径类型
@pytest.mark.parametrize("path_type", [str, Path])
def test_path_types(tmp_path, path_type):
    """测试接受字符串和Path对象作为输入"""
    test_dir = tmp_path / "path_type_test"
    result = create_swanlog_dir(path_type(test_dir))

    assert result == path_type(test_dir)
    assert os.path.exists(test_dir)


# 测试空目录创建 .gitignore 但目录非空的情况
def test_gitignore_in_non_empty_directory(tmp_path):
    """测试在非空目录中不创建 .gitignore"""
    test_dir = tmp_path / "non_empty_swanlog"
    test_dir.mkdir()

    # 创建一些文件
    (test_dir / "file1.txt").touch()
    (test_dir / "file2.txt").touch()

    create_swanlog_dir(test_dir)

    gitignore = test_dir / ".gitignore"
    assert not gitignore.exists()  # 不应创建 .gitignore


# 测试 .gitignore 文件内容
def test_gitignore_content(tmp_path):
    """测试 .gitignore 文件内容是否正确"""
    test_dir = tmp_path / "gitignore_test"
    test_dir.mkdir()

    create_swanlog_dir(test_dir)

    gitignore = test_dir / ".gitignore"
    assert gitignore.exists()
    with open(gitignore, "r", encoding="utf-8") as f:
        content = f.read()
        assert content == "*"
        assert len(content) == 1  # 只有一个字符 '*'


# 测试在已存在 .gitignore 的情况下不覆盖
def test_existing_gitignore(tmp_path):
    """测试目录中已存在 .gitignore 时不覆盖"""
    test_dir = tmp_path / "existing_gitignore"
    test_dir.mkdir()

    # 预先创建 .gitignore 文件
    gitignore = test_dir / ".gitignore"
    with open(gitignore, "w", encoding="utf-8") as f:
        f.write("custom_content")

    create_swanlog_dir(test_dir)

    # 验证文件内容未被覆盖
    with open(gitignore, "r", encoding="utf-8") as f:
        assert f.read() == "custom_content"
