#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/7/18 15:28
@File: test_env.py
@IDE: pycharm
@Description:
    测试swanlab.env模块
"""
import os
from pathlib import Path

import pytest

import swanlab
from swanlab.env import SwanLabEnv, is_interactive, get_save_dir, remove_host_suffix, create_swanlog_dir


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
    def mock_makedirs(path, exist_ok=False):
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
