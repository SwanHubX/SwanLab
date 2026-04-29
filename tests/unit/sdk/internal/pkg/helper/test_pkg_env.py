"""
@author: cunyue
@file: test_pkg_env.py
@time: 2026/4/29
@description: 测试 is_interactive / is_jupyter / is_jupyter_supports_stdin
"""

from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.pkg.helper.env import is_interactive, is_jupyter_supports_stdin


class TestJupyterSupportsStdin:
    """测试 is_jupyter_supports_stdin 函数"""

    def test_kernel_allows_stdin(self):
        """kernel._allow_stdin=True 时返回 True"""
        ip = MagicMock()
        ip.kernel._allow_stdin = True
        with patch("IPython.core.getipython.get_ipython", return_value=ip):
            assert is_jupyter_supports_stdin() is True

    def test_kernel_disallows_stdin(self):
        """kernel._allow_stdin=False 时返回 False（如 nbconvert --execute）"""
        ip = MagicMock()
        ip.kernel._allow_stdin = False
        with patch("IPython.core.getipython.get_ipython", return_value=ip):
            assert is_jupyter_supports_stdin() is False

    def test_no_kernel_attribute(self):
        """没有 kernel 属性时保守返回 True"""
        ip = MagicMock(spec=[])  # 空 spec，无任何属性
        with patch("IPython.core.getipython.get_ipython", return_value=ip):
            assert is_jupyter_supports_stdin() is True

    def test_exception_returns_true(self):
        """异常时保守返回 True"""
        with patch("IPython.core.getipython.get_ipython", side_effect=RuntimeError):
            assert is_jupyter_supports_stdin() is True


class TestIsInteractiveJupyterStdin:
    """测试 is_interactive 对 Jupyter stdin 的处理"""

    def test_jupyter_nbconvert_not_interactive(self, monkeypatch):
        """nbconvert 环境下 is_interactive 应返回 False"""
        monkeypatch.setattr("swanlab.sdk.internal.pkg.helper.env.is_jupyter", lambda: True)
        monkeypatch.setattr("swanlab.sdk.internal.pkg.helper.env.is_jupyter_supports_stdin", lambda: False)
        # mock stdin 有 fileno 且非 tty
        mock_stdin = MagicMock()
        mock_stdin.fileno.return_value = 0
        monkeypatch.setattr("sys.stdin", mock_stdin)
        monkeypatch.setattr("os.isatty", lambda fd: False)
        assert is_interactive() is False

    def test_jupyter_interactive_notebook(self, monkeypatch):
        """正常 Jupyter Notebook 下 is_interactive 应返回 True"""
        monkeypatch.setattr("swanlab.sdk.internal.pkg.helper.env.is_jupyter", lambda: True)
        monkeypatch.setattr("swanlab.sdk.internal.pkg.helper.env.is_jupyter_supports_stdin", lambda: True)
        mock_stdin = MagicMock()
        mock_stdin.fileno.return_value = 0
        monkeypatch.setattr("sys.stdin", mock_stdin)
        monkeypatch.setattr("os.isatty", lambda fd: False)
        assert is_interactive() is True

    def test_tty_always_interactive(self, monkeypatch):
        """TTY 环境始终为可交互"""
        mock_stdin = MagicMock()
        mock_stdin.fileno.return_value = 0
        monkeypatch.setattr("sys.stdin", mock_stdin)
        monkeypatch.setattr("os.isatty", lambda fd: True)
        assert is_interactive() is True
