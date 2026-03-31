"""
@author: cunyue
@file: test_runtime.py
@time: 2026/3/31
@description: Tests for runtime information collection
"""

import sys
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.run.system.environment import runtime
from swanlab.sdk.typings.run.system import RuntimeSnapshot


class TestGetOs:
    """Tests for runtime.get_os()"""

    @pytest.mark.parametrize(
        "platform_return,expected",
        [
            ("macOS-14.0-arm64-arm-64bit", "macOS-14.0-arm64-arm-64bit"),
            ("Linux-5.15.0-x86_64-with-glibc2.35", "Linux-5.15.0-x86_64-with-glibc2.35"),
        ],
    )
    def test_get_os_success(self, platform_return, expected):
        """成功获取操作系统平台信息"""
        with patch("platform.platform", return_value=platform_return):
            assert runtime.get_os() == expected

    def test_get_os_returns_none_on_error(self):
        """platform.platform 异常时返回 None"""
        with patch("platform.platform", side_effect=RuntimeError("mock error")):
            assert runtime.get_os() is None


class TestGetOsPretty:
    """Tests for runtime.get_os_pretty()"""

    @pytest.mark.parametrize(
        "os_release,expected",
        [
            ({"PRETTY_NAME": "Ubuntu 22.04.3 LTS"}, "Ubuntu 22.04.3 LTS"),
            ({"NAME": "Ubuntu"}, None),
        ],
    )
    def test_get_os_pretty(self, os_release, expected):
        """测试获取操作系统友好名称"""
        with patch("platform.freedesktop_os_release", return_value=os_release):
            assert runtime.get_os_pretty() == expected

    def test_get_os_pretty_returns_none_on_error(self):
        """freedesktop_os_release 异常时返回 None"""
        with patch("platform.freedesktop_os_release", side_effect=OSError("not available")):
            assert runtime.get_os_pretty() is None


class TestGetHostname:
    """Tests for runtime.get_hostname()"""

    def test_get_hostname_success(self):
        """成功获取主机名"""
        with patch("socket.gethostname", return_value="my-macbook-pro"):
            assert runtime.get_hostname() == "my-macbook-pro"

    def test_get_hostname_returns_none_on_error(self):
        """socket.gethostname 异常时返回 None"""
        with patch("socket.gethostname", side_effect=OSError("network error")):
            assert runtime.get_hostname() is None


class TestGetPid:
    """Tests for runtime.get_pid()"""

    def test_get_pid_success(self):
        """成功获取进程 ID"""
        with patch("os.getpid", return_value=12345):
            assert runtime.get_pid() == 12345

    def test_get_pid_returns_none_on_error(self):
        """os.getpid 异常时返回 None"""
        with patch("os.getpid", side_effect=RuntimeError("mock error")):
            assert runtime.get_pid() is None


class TestGetCwd:
    """Tests for runtime.get_cwd()"""

    def test_get_cwd_success(self):
        """成功获取当前工作目录"""
        with patch("os.getcwd", return_value="/home/user/project"):
            assert runtime.get_cwd() == "/home/user/project"

    def test_get_cwd_returns_none_on_error(self):
        """os.getcwd 异常时返回 None"""
        with patch("os.getcwd", side_effect=OSError("directory removed")):
            assert runtime.get_cwd() is None


class TestGetPythonVersion:
    """Tests for runtime.get_python_version()"""

    def test_get_python_version_success(self):
        """成功获取 Python 版本"""
        with patch("platform.python_version", return_value="3.11.5"):
            assert runtime.get_python_version() == "3.11.5"

    def test_get_python_version_returns_none_on_error(self):
        """platform.python_version 异常时返回 None"""
        with patch("platform.python_version", side_effect=RuntimeError("mock error")):
            assert runtime.get_python_version() is None


class TestGetCommand:
    """Tests for runtime.get_command()"""

    @pytest.mark.parametrize(
        "system,argv,expected",
        [
            ("Darwin", ["python", "script.py", "--arg1", "value1"], "python script.py --arg1 value1"),
            ("Windows", ["python.exe", "main.py"], "python.exe main.py"),
        ],
    )
    def test_get_command_non_linux(self, system, argv, expected):
        """非 Linux 系统使用 sys.argv"""
        with (
            patch("platform.system", return_value=system),
            patch.object(sys, "argv", argv),
        ):
            assert runtime.get_command() == expected

    def test_get_command_linux_success(self):
        """Linux 系统从 /proc/self/cmdline 读取"""
        cmdline_content = b"python\x00script.py\x00--arg\x00value\x00"
        mock_file = MagicMock()
        mock_file.read.return_value = cmdline_content
        mock_file.__enter__ = MagicMock(return_value=mock_file)
        mock_file.__exit__ = MagicMock(return_value=False)

        with (
            patch("platform.system", return_value="Linux"),
            patch("builtins.open", return_value=mock_file),
        ):
            assert runtime.get_command() == "python script.py --arg value"

    @pytest.mark.parametrize(
        "linux_result,argv,expected",
        [
            (None, ["python", "fallback.py"], "python fallback.py"),
            ("", ["python", "fallback.py"], "python fallback.py"),
        ],
    )
    def test_get_command_linux_fallback(self, linux_result, argv, expected):
        """Linux 系统 cmdline 读取失败时回退到 sys.argv"""
        with (
            patch("platform.system", return_value="Linux"),
            patch("swanlab.sdk.internal.run.system.environment.runtime._get_command_linux", return_value=linux_result),
            patch.object(sys, "argv", argv),
        ):
            assert runtime.get_command() == expected


class TestGetCommandLinux:
    """Tests for runtime._get_command_linux()"""

    @pytest.mark.parametrize(
        "cmdline,expected",
        [
            (b"/usr/bin/python3\x00/home/user/app.py\x00--debug\x00", "/usr/bin/python3 /home/user/app.py --debug"),
            (b"python\x00script.py\x00--message\x00hello world\x00", "python script.py --message hello world"),
            ("python 脚本.py --中文 参数\x00".encode("utf-8"), "python 脚本.py --中文 参数"),
        ],
    )
    def test_get_command_linux_success(self, cmdline, expected):
        """正常读取 Linux cmdline"""
        mock_file = MagicMock()
        mock_file.read.return_value = cmdline
        mock_file.__enter__ = MagicMock(return_value=mock_file)
        mock_file.__exit__ = MagicMock(return_value=False)

        with patch("builtins.open", return_value=mock_file):
            assert runtime._get_command_linux() == expected

    @pytest.mark.parametrize(
        "side_effect",
        [
            FileNotFoundError("proc not mounted"),
            OSError("permission denied"),
        ],
    )
    def test_get_command_linux_error(self, side_effect):
        """读取 /proc/self/cmdline 失败时返回 None"""
        with patch("builtins.open", side_effect=side_effect):
            assert runtime._get_command_linux() is None


class TestRuntimeGet:
    """Tests for runtime.get()"""

    def test_get_returns_runtime_snapshot(self):
        """get() 返回 RuntimeSnapshot 实例"""
        with (
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_os", return_value="macOS-14.0"),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_os_pretty", return_value="macOS Sonoma"),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_hostname", return_value="my-mac"),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_pid", return_value=12345),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_cwd", return_value="/home/user"),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_python_version", return_value="3.11.0"),
            patch(
                "swanlab.sdk.internal.run.system.environment.runtime.get_python_verbose", return_value="3.11.0 (main)"
            ),
            patch(
                "swanlab.sdk.internal.run.system.environment.runtime.get_python_executable",
                return_value="/usr/bin/python3",
            ),
            patch("swanlab.sdk.internal.run.system.environment.runtime.get_command", return_value="python app.py"),
        ):
            result = runtime.get()

        assert isinstance(result, RuntimeSnapshot)
        assert result.os == "macOS-14.0"
        assert result.hostname == "my-mac"
        assert result.pid == 12345

    def test_get_with_none_values(self):
        """部分信息获取失败时仍返回 RuntimeSnapshot"""
        mocks = [
            "get_os",
            "get_os_pretty",
            "get_hostname",
            "get_pid",
            "get_cwd",
            "get_python_version",
            "get_python_verbose",
            "get_python_executable",
            "get_command",
        ]
        with patch("swanlab.sdk.internal.run.system.environment.runtime") as mock_runtime:
            for m in mocks:
                patcher = patch.object(mock_runtime, m, return_value=None)
                patcher.start()

            # 直接调用 get 函数测试 None 值处理
            result = runtime.get()

        assert isinstance(result, RuntimeSnapshot)
