"""
@author: cunyue
@file: test_toolkit_logger.py
@time: 2025/10/10 22:27
@description: 测试工具中的日志模块
"""

import nanoid

from swanlab.toolkit import SwanKitLogger


class TestSwanKitLog:
    """
    测试日志模块
    """

    def test_invalid_level(self):
        """
        测试无效日志等级，默认设置为 info
        """
        name = nanoid.generate()
        t = SwanKitLogger(name, level="invalid_level")
        assert t.level == "info"

    def test_enable_default(self, capsys):
        """
        测试开启日志
        """
        levels = ["debug", "info", "warning", "error", "critical"]
        for level in levels:
            name = nanoid.generate()
            text = nanoid.generate()
            t = SwanKitLogger(name, level=level)
            for le in levels:
                getattr(t, le)(text)
                out, err = capsys.readouterr()
                if levels.index(le) >= levels.index(level):
                    assert text in out
                    assert name in out
                    assert err == ""
                else:
                    assert out == ""
                    assert err == ""

    def test_disable(self, capsys):
        """
        测试关闭日志
        """
        levels = ["debug", "info", "warning", "error", "critical"]
        for level in levels:
            name = nanoid.generate()
            text = nanoid.generate()
            t = SwanKitLogger(name, level=level)
            t.disable_log()
            for le in levels:
                getattr(t, le)(text)
                out, err = capsys.readouterr()
                assert out == ""
                assert err == ""

    def test_set_level(self, capsys):
        """
        测试设置日志等级
        """
        levels = ["debug", "info", "warning", "error", "critical"]
        for le in levels:
            name = nanoid.generate()
            text = nanoid.generate()
            t = SwanKitLogger(name, level="debug")
            t.level = le
            getattr(t, le)(text)
            out, err = capsys.readouterr()
            assert text in out
            assert name in out
            assert err == ""

    def test_markup_disabled_by_default(self, capsys):
        """
        测试默认不解析用户日志中的 Rich markup
        """
        name = nanoid.generate()
        text = "[/tmp/e2e/sub/compose.yml]"
        t = SwanKitLogger(name)

        t.info(text)
        out, err = capsys.readouterr()

        assert text in out
        assert name in out
        assert err == ""

    def test_markup_kwarg_can_be_explicitly_passed(self, capsys):
        """
        测试显式传递 markup 参数不会产生重复参数错误
        """
        name = nanoid.generate()
        t = SwanKitLogger(name)

        t.info("[bold]hello[/bold]", markup=True)
        out, err = capsys.readouterr()

        assert "hello" in out
        assert err == ""
