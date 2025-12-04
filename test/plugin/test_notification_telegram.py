#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 TelegramCallback 插件
参考 Discord/Slack 的测试风格，使用 responses 库 mock HTTP 请求
"""
import pytest
import responses
from unittest.mock import patch, MagicMock
from swanlab.plugin import TelegramCallback
from swanlab.plugin.notification import TelegramBot


class TestTelegramBot:
    """测试 TelegramBot 类"""

    def test_init(self):
        """测试初始化"""
        bot = TelegramBot("test_token", "123456")
        assert bot.bot_token == "test_token"
        assert bot.chat_id == "123456"
        assert bot.api_base == "https://api.telegram.org/bottest_token"

    @responses.activate
    def test_send_message_success(self):
        """测试发送消息成功（使用 responses mock）"""
        # Mock Telegram API 响应
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 123}},
            status=200,
        )

        bot = TelegramBot("test_token", "123456")
        result = bot.send_message("Hello, World!")

        assert result["ok"] is True
        assert result["result"]["message_id"] == 123
        assert len(responses.calls) == 1
        
        # 验证请求参数
        import json
        request_body = json.loads(responses.calls[0].request.body)
        assert request_body["chat_id"] == "123456"
        assert request_body["text"] == "Hello, World!"
        assert request_body["parse_mode"] == "HTML"

    @responses.activate
    def test_send_message_with_markdown(self):
        """测试使用 Markdown 格式发送消息"""
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 124}},
            status=200,
        )

        bot = TelegramBot("test_token", "123456")
        result = bot.send_message("**Bold** text", parse_mode="Markdown")

        assert result["ok"] is True
        import json
        request_body = json.loads(responses.calls[0].request.body)
        assert request_body["parse_mode"] == "Markdown"

    @responses.activate
    def test_send_message_api_error(self):
        """测试 API 返回错误"""
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": False, "description": "Bad Request: chat not found"},
            status=200,
        )

        bot = TelegramBot("test_token", "invalid_chat")
        result = bot.send_message("Hello")

        assert result["ok"] is False
        assert "chat not found" in result["description"]

    @responses.activate
    def test_send_message_network_error(self):
        """测试网络错误"""
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            body=Exception("Connection error"),
        )

        bot = TelegramBot("test_token", "123456")
        with pytest.raises(Exception):
            bot.send_message("Hello")


class TestTelegramCallback:
    """测试 TelegramCallback 类"""

    def test_init(self):
        """测试初始化"""
        callback = TelegramCallback("test_token", "123456", language="zh")
        assert callback.bot.bot_token == "test_token"
        assert callback.bot.chat_id == "123456"
        assert callback.language == "zh"
        assert callback.notify_on_start is False

    def test_init_with_notify_on_start(self):
        """测试带 notify_on_start 参数的初始化"""
        callback = TelegramCallback("test_token", "123456", notify_on_start=True)
        assert callback.notify_on_start is True

    def test_init_english(self):
        """测试英文语言初始化"""
        callback = TelegramCallback("test_token", "123456", language="en")
        assert callback.language == "en"

    def test_str(self):
        """测试 __str__ 方法"""
        callback = TelegramCallback("test_token", "123456")
        assert str(callback) == "TelegramCallback"

    def test_on_init(self):
        """测试 on_init 方法"""
        callback = TelegramCallback("test_token", "123456")
        callback.on_init("test_project", "test_workspace")
        assert callback.project == "test_project"
        assert callback.workspace == "test_workspace"

    def test_on_init_with_optional_params(self):
        """测试 on_init 方法带可选参数"""
        callback = TelegramCallback("test_token", "123456")
        callback.on_init("test_project", "test_workspace", public=True, logdir="/tmp/logs")
        assert callback.project == "test_project"
        assert callback.workspace == "test_workspace"

    def test_before_init_experiment(self):
        """测试 before_init_experiment 方法"""
        callback = TelegramCallback("test_token", "123456")
        callback.before_init_experiment(
            run_id="test_run_id",
            exp_name="test_exp",
            description="test description",
            colors=("#fff", "#000"),
        )
        assert callback.run_id == "test_run_id"
        assert callback.exp_name == "test_exp"
        assert callback.description == "test description"

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_success_zh(self, mock_get_url):
        """测试创建成功通知内容（中文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error=None)

        assert "SwanLab 消息通知" in content
        assert "实验已成功完成" in content
        assert "test_project" in content
        assert "test_workspace" in content
        assert "test_exp" in content
        assert "https://swanlab.cn/exp/123" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_success_en(self, mock_get_url):
        """测试创建成功通知内容（英文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="en")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error=None)

        assert "SwanLab Notification" in content
        assert "Experiment completed successfully" in content
        assert "Project:" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_error_zh(self, mock_get_url):
        """测试创建错误通知内容（中文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error="内存不足")

        assert "实验遇到错误" in content
        assert "内存不足" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_error_en(self, mock_get_url):
        """测试创建错误通知内容（英文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="en")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error="Out of memory")

        assert "Experiment failed" in content
        assert "Out of memory" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_offline_zh(self, mock_get_url):
        """测试离线模式通知内容（中文）"""
        mock_get_url.return_value = None

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error=None)

        assert "离线模式" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_offline_en(self, mock_get_url):
        """测试离线模式通知内容（英文）"""
        mock_get_url.return_value = None

        callback = TelegramCallback("test_token", "123456", language="en")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="stop", error=None)

        assert "offline mode" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_start_zh(self, mock_get_url):
        """测试创建开始通知内容（中文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="start")

        assert "实验已开始" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_start_en(self, mock_get_url):
        """测试创建开始通知内容（英文）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="en")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        content = callback._create_content(event="start")

        assert "Experiment started" in content

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_create_content_none_description(self, mock_get_url):
        """测试 description 为 None 时的处理"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = None

        content = callback._create_content(event="stop", error=None)

        assert "N/A" in content  # None 应该被替换为 N/A

    @responses.activate
    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_send_msg_success(self, mock_get_url):
        """测试 send_msg 成功发送"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 123}},
            status=200,
        )

        callback = TelegramCallback("test_token", "123456")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        # 不应该抛出异常
        callback.send_msg("Test message")
        assert len(responses.calls) == 1

    @responses.activate
    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_send_msg_failure(self, mock_get_url):
        """测试 send_msg 发送失败"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": False, "description": "Bad Request"},
            status=200,
        )

        callback = TelegramCallback("test_token", "123456")
        # 不应该抛出异常，只打印错误信息
        callback.send_msg("Test message")
        assert len(responses.calls) == 1

    @responses.activate
    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_on_stop_integration(self, mock_get_url):
        """测试 on_stop 完整流程（集成测试）"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 123}},
            status=200,
        )

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        callback.on_stop(error=None)

        # 验证发送了请求
        assert len(responses.calls) == 1
        import json
        request_body = json.loads(responses.calls[0].request.body)
        assert "实验已成功完成" in request_body["text"]

    @responses.activate
    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_on_stop_with_error_integration(self, mock_get_url):
        """测试 on_stop 带错误的完整流程"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 123}},
            status=200,
        )

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        callback.on_stop(error="CUDA out of memory")

        assert len(responses.calls) == 1
        import json
        request_body = json.loads(responses.calls[0].request.body)
        assert "实验遇到错误" in request_body["text"]
        assert "CUDA out of memory" in request_body["text"]

    @responses.activate
    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_on_run_with_notify_integration(self, mock_get_url):
        """测试 on_run 带通知的完整流程"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"
        responses.add(
            responses.POST,
            "https://api.telegram.org/bottest_token/sendMessage",
            json={"ok": True, "result": {"message_id": 123}},
            status=200,
        )

        callback = TelegramCallback("test_token", "123456", notify_on_start=True, language="zh")
        callback.project = "test_project"
        callback.workspace = "test_workspace"
        callback.exp_name = "test_exp"
        callback.description = "test description"

        callback.on_run()

        assert len(responses.calls) == 1
        import json
        request_body = json.loads(responses.calls[0].request.body)
        assert "实验已开始" in request_body["text"]

    @responses.activate
    def test_on_run_without_notify(self):
        """测试 on_run 方法（禁用开始通知）"""
        callback = TelegramCallback("test_token", "123456", notify_on_start=False)

        callback.on_run()

        # 不应该发送任何请求
        assert len(responses.calls) == 0


class TestTelegramCallbackEdgeCases:
    """测试边界情况"""

    def test_empty_bot_token(self):
        """测试空 bot_token"""
        callback = TelegramCallback("", "123456")
        assert callback.bot.bot_token == ""
        assert callback.bot.api_base == "https://api.telegram.org/bot"

    def test_special_characters_in_chat_id(self):
        """测试 chat_id 包含特殊字符（如频道 @username）"""
        callback = TelegramCallback("test_token", "@my_channel")
        assert callback.bot.chat_id == "@my_channel"

    @patch("swanlab.plugin.notification.swanlab.get_url")
    def test_html_escape_in_content(self, mock_get_url):
        """测试内容中的 HTML 特殊字符"""
        mock_get_url.return_value = "https://swanlab.cn/exp/123"

        callback = TelegramCallback("test_token", "123456", language="zh")
        callback.project = "test<project>"
        callback.workspace = "test&workspace"
        callback.exp_name = "test\"exp\""
        callback.description = "test'description'"

        # 不应该抛出异常
        content = callback._create_content(event="stop", error=None)
        assert "test<project>" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

