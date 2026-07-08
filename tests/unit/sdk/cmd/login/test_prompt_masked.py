"""
@author: opencode
@file: test_prompt_masked.py
@time: 2026/7/8
@description: prompt_masked 星号遮罩输入辅助函数单测（基于 pwinput）
"""

import sys
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.cmd.utils import prompt_masked


class TestPromptMasked:
    """通过 mock pwinput.getch 验证遮罩输入与 Ctrl+C/Ctrl+D 拦截"""

    @staticmethod
    def _setup(monkeypatch, chars: str):
        """让 pwinput 走遮罩分支（sys.stdin is sys.__stdin__），并注入逐字符输入。"""
        import pwinput

        # 令 pwinput 不回退到 getpass：sys.stdin 必须是 sys.__stdin__
        monkeypatch.setattr("sys.stdin", sys.__stdin__)
        it = iter(chars)

        def _fake_getch():
            return next(it)

        monkeypatch.setattr(pwinput, "getch", _fake_getch)
        # 捕获 stdout 输出
        output: list[str] = []
        monkeypatch.setattr("sys.stdout", MagicMock(write=output.append, flush=lambda: None))
        return output

    def test_basic_input(self, monkeypatch):
        """逐字符输入 abc + 回车，返回 'abc'，屏幕输出遮罩符号 + 换行"""
        output = self._setup(monkeypatch, "abc\r")
        result = prompt_masked()
        assert result == "abc"
        out = "".join(output)
        assert "***\n" in out  # 遮罩 + 回车
        assert "\x1b[?2004l" in out  # 禁用括号化粘贴
        assert "\x1b[?2004h" in out  # 恢复括号化粘贴

    def test_backspace(self, monkeypatch):
        """输入 ab + 退格 + c + 回车 → 返回 'ac'"""
        self._setup(monkeypatch, "ab\x7fc\r")
        result = prompt_masked()
        assert result == "ac"

    def test_backspace_empty_buffer(self, monkeypatch):
        """空缓冲区下退格不应报错"""
        self._setup(monkeypatch, "\x7f\x7fa\r")
        result = prompt_masked()
        assert result == "a"

    def test_ctrl_c_raises(self, monkeypatch):
        """Ctrl+C 抛 KeyboardInterrupt（pwinput 原生会吞掉，wrapper 必须拦截）"""
        self._setup(monkeypatch, "\x03")
        with pytest.raises(KeyboardInterrupt):
            prompt_masked()

    def test_ctrl_d_raises_eof(self, monkeypatch):
        """Ctrl+D 抛 EOFError"""
        self._setup(monkeypatch, "\x04")
        with pytest.raises(EOFError):
            prompt_masked()

    def test_ctrl_z_raises_eof(self, monkeypatch):
        """Ctrl+Z (Windows EOF) 抛 EOFError"""
        self._setup(monkeypatch, "\x1a")
        with pytest.raises(EOFError):
            prompt_masked()

    def test_custom_mask(self, monkeypatch):
        """自定义遮罩符号"""
        output = self._setup(monkeypatch, "ab\r")
        result = prompt_masked(mask="•")
        assert result == "ab"
        out = "".join(output)
        assert "••\n" in out


class TestPromptMasked_Fallback:
    """pwinput 在非 tty 环境（sys.stdin is not sys.__stdin__）回退到 getpass"""

    def test_fallback_to_getpass(self, monkeypatch):
        """sys.stdin 被替换时 pwinput 回退 getpass，Ctrl+C/Ctrl+D 由 getpass 原生处理"""
        import getpass

        monkeypatch.setattr("sys.stdin", MagicMock())  # 不是 sys.__stdin__
        called = {"value": False}

        def _fake_getpass(prompt=""):
            called["value"] = True
            return "fallback-key"

        monkeypatch.setattr(getpass, "getpass", _fake_getpass)
        monkeypatch.setattr("sys.stdout", MagicMock(write=lambda _: None, flush=lambda: None))
        result = prompt_masked()
        assert called["value"] is True
        assert result == "fallback-key"
