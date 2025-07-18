#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/21 16:57
@File: pytest_log.py
@IDE: pycharm
@Description:
    测试swanlog类，只需测试其日志监听功能
"""
import os
import sys
import time

import pytest
from nanoid import generate

from swanlab.log import swanlog
from swanlab.log.log import clean_control_chars, remove_control_sequences
from swanlab.log.type import LogData, ProxyType
from tutils import TEMP_PATH


class TestSwanLogInstall:
    """
    目前在设计上不希望外界实例化SwanLog，所以不提供实例化测试
    """

    @staticmethod
    def teardown_method():
        # 每个测试方法后执行
        try:
            swanlog.reset()
        except RuntimeError:
            pass

    @staticmethod
    def setup_method():
        # 每个测试方法前执行
        try:
            swanlog.reset()
        except RuntimeError:
            pass

    @staticmethod
    def start_proxy(proxy_type: ProxyType = "all", max_log_length=1024):
        """
        新建一个回调函数
        """
        console_dir = os.path.join(TEMP_PATH, str(generate()))
        os.mkdir(console_dir)
        # 日志文件路径
        from datetime import datetime

        now = datetime.now().strftime("%Y-%m-%d")
        log_file = os.path.join(console_dir, f"{now}.log")

        # 拿到日志文件句柄
        def write_handler(log_data: LogData):
            with open(log_file, "a") as f:
                for content in log_data["contents"]:
                    f.write(content["message"] + "\n")

        swanlog.start_proxy(proxy_type, max_log_length, write_handler)
        if proxy_type == 'all':
            assert getattr(swanlog, '_SwanLog__origin_stdout_write') is not None
            assert getattr(swanlog, '_SwanLog__origin_stderr_write') is not None
        elif proxy_type == 'stdout':
            assert getattr(swanlog, '_SwanLog__origin_stdout_write') is not None
            assert getattr(swanlog, '_SwanLog__origin_stderr_write') is None
        elif proxy_type == 'stderr':
            assert getattr(swanlog, '_SwanLog__origin_stdout_write') is None
            assert getattr(swanlog, '_SwanLog__origin_stderr_write') is not None
        else:
            raise ValueError(f"Invalid proxy type: {proxy_type}")
        return log_file

    def test_global_install(self):
        self.start_proxy()
        assert swanlog.proxied is True
        with pytest.raises(RuntimeError) as e:
            self.start_proxy()
        assert str(e.value) == "Std Proxy is already started"
        swanlog.stop_proxy()
        assert swanlog.proxied is False
        self.start_proxy()
        assert swanlog.proxied is True

    def test_write_after_uninstall(self):
        """
        在卸载后打印，此时应该不会写入日志文件
        """
        log_file = self.start_proxy()
        swanlog.stop_proxy()
        print("\ntest write after uninstall")
        a = generate()
        print(a)
        b = generate()
        print(b)
        assert not os.path.exists(log_file)

    def test_write_to_file(self):
        """
        测试写入日志到文件
        """
        log_file = self.start_proxy()
        print("test write to file")
        a = generate()
        print(a)
        b = generate()
        print(b)
        assert os.path.exists(log_file)
        # 比较最后两行内容
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[-2] == a + "\n"
            assert content[-1] == b + "\n"
        # 卸载后再次 print，此时应该不会写入日志文件
        swanlog.stop_proxy()
        print("\ntest after stop proxy")
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[-2] == a + "\n"
            assert content[-1] == b + "\n"

    def test_write_to_file_long_test(self):
        log_file = self.start_proxy()
        # 获取默认最大长度
        max_len = 1024
        # 默认最大长度为1024
        a = generate(size=3000)
        print(a)
        time.sleep(0.1)
        with open(os.path.join(log_file), "r") as f:
            content = f.readlines()
            assert content[-1] == a[:max_len] + "\n"

    def test_write_logging_to_file(self):
        # FIXME 不知道为什么此函数在 pycharm 的测试中如果不设置路径为 ./test/unit 而是 ./test 会报错
        # FIXME 在云端测试模式、联合测试环境下，硬件监控的定时器好像没有在其他地方取消，这会导致设置 debug 级别的同时捕获硬件监控的日志导致报错
        log_file = self.start_proxy()
        swanlog.level = 'warning'
        print("test write to file")
        a = generate()
        swanlog.warning(a)
        b = generate()
        swanlog.error(b)
        time.sleep(0.1)
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[-2] == "swanlab: " + a + "\n"
            assert content[-1] == "swanlab: " + b + "\n"

    def test_can_write_logging(self):
        # FIXME 不知道为什么此函数在 pycharm 的测试中如果不设置路径为 ./test/unit 而是 ./test 会报错
        log_file = self.start_proxy()
        print("test write to file")
        a = generate()
        # debug 默认不打印
        swanlog.debug(a)
        b = generate()
        swanlog.info(b)
        time.sleep(0.1)
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[-2] == "test write to file\n"
            assert content[-1] == "swanlab: " + b + "\n"

    def test_stderr_write(self):
        """
        测试stderr的写入
        """
        log_file = self.start_proxy()
        print("test write to file")
        a = generate()
        sys.stderr.write(a + "\n")
        b = generate()
        sys.stderr.write(b + "\n")
        assert os.path.exists(log_file)
        time.sleep(0.1)
        # 比较最后两行内容
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[0] == "test write to file\n"
            assert content[1] == a + "\n"
            assert content[2] == b + "\n"

    def test_stderr_only(self):
        log_file = self.start_proxy('stderr')
        print("test write to file")
        a = generate()
        sys.stderr.write(a + "\n")
        b = generate()
        sys.stderr.write(b + "\n")
        time.sleep(0.1)
        assert os.path.exists(log_file)
        # 比较最后两行内容
        with open(log_file, "r") as f:
            content = f.readlines()
            assert content[0] == a + "\n"
            assert content[1] == b + "\n"

    # def test_write_sharding(self, monkeypatch):
    #     """
    #     测试日志文件分片
    #     """
    #     console_dir = self.create_console_dir()
    #     with freeze_time('2020-10-06'):
    #         swanlog.install(console_dir)
    #         print("1234")
    #         assert os.path.exists(os.path.join(console_dir, "2020-10-06.log"))
    #     with freeze_time('2020-10-07'):
    #         p = os.path.join(console_dir, "2020-10-07.log")
    #         assert not os.path.exists(p)
    #         print("1234")
    #         assert os.path.exists(p)


class TestCleanControlChars:
    @staticmethod
    def test_handles_text_with_no_newline():
        text = "Hello World"
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == []
        assert buffer == "Hello World"

    @staticmethod
    def test_removes_control_sequences_and_returns_cleaned_lines():
        text = "Hello\rWorld\nLine\x1b[AOne\nLine Two"
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == ["World", "One"]
        assert buffer == "Line Two"

    @staticmethod
    def test_handles_empty_text_and_returns_empty():
        text = ""
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == []
        assert buffer == ""

    @staticmethod
    def test_handles_text_with_only_newline():
        text = "\n"
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == []
        assert buffer == ""

    @staticmethod
    def test_handles_text_with_multiple_newlines():
        text = "Line One\nLine Two\nLine Three\n"
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == ["Line One", "Line Two", "Line Three"]
        assert buffer == ""

    @staticmethod
    def test_removes_ansi_escape_sequences():
        text = "Line\x1b[31mRed\x1b[0mOne\nLine\x1b[32mGreen\x1b[0mTwo\n"
        cleaned_lines, buffer = clean_control_chars(text)
        assert cleaned_lines == ["LineRedOne", "LineGreenTwo"]
        assert buffer == ""


class TestRemoveControlSequences:
    @staticmethod
    def test_removes_carriage_return_and_returns_content_after():
        line = "Hello\rWorld"
        assert remove_control_sequences(line) == "World"

    @staticmethod
    def removes_escape_sequence_and_returns_content_after():
        line = "Hello\x1b[AWorld"
        assert remove_control_sequences(line) == "World"

    @staticmethod
    def returns_original_line_if_no_control_sequence():
        line = "Hello World"
        assert remove_control_sequences(line) == "Hello World"

    @staticmethod
    def handles_multiple_control_sequences_and_returns_content_after_last():
        line = "Hello\rWorld\x1b[AUniverse"
        assert remove_control_sequences(line) == "Universe"

    @staticmethod
    def handles_empty_string_and_returns_empty():
        line = ""
        assert remove_control_sequences(line) == ""

    @staticmethod
    def handles_only_control_sequence_and_returns_empty():
        line = "\r"
        assert remove_control_sequences(line) == ""
