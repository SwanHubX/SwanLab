"""
@author: cunyue
@file: test_log.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.cmd.log 中的函数
"""

from unittest.mock import MagicMock

from swanlab.sdk.cmd.log import log, log_text


class TestLog:
    def test_log_calls_run_log(self, monkeypatch):
        """log() 应调用 run.log() 并传递参数"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.helper.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.log.get_run", lambda: mock_run)

        log({"loss": 0.5, "accuracy": 0.95})

        mock_run.log.assert_called_once_with({"loss": 0.5, "accuracy": 0.95}, None)

    def test_log_with_step(self, monkeypatch):
        """log() 应正确传递 step 参数"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.helper.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.log.get_run", lambda: mock_run)

        log({"loss": 0.5}, step=10)

        mock_run.log.assert_called_once_with({"loss": 0.5}, 10)


class TestLogText:
    def test_log_text_calls_run_log_text(self, monkeypatch):
        """log_text() 应调用 run.log_text() 并传递参数"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.helper.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.log.get_run", lambda: mock_run)

        log_text("output", "Training started")

        mock_run.log_text.assert_called_once_with("output", "Training started", None, None)

    def test_log_text_with_caption(self, monkeypatch):
        """log_text() 应正确传递 caption 参数"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.helper.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.log.get_run", lambda: mock_run)

        log_text("prediction", "cat", caption="Model output")

        mock_run.log_text.assert_called_once_with("prediction", "cat", "Model output", None)

    def test_log_text_with_step(self, monkeypatch):
        """log_text() 应正确传递 step 参数"""
        mock_run = MagicMock()
        monkeypatch.setattr("swanlab.sdk.cmd.helper.has_run", lambda: True)
        monkeypatch.setattr("swanlab.sdk.cmd.log.get_run", lambda: mock_run)

        log_text("status", "Epoch complete", step=5)

        mock_run.log_text.assert_called_once_with("status", "Epoch complete", None, 5)
