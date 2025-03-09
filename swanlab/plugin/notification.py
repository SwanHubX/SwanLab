"""
Notification plugin for SwanLab.
"""

from swankit.callback import SwanKitCallback
import swanlab
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Optional, Dict, Any
import hmac
import hashlib
import base64
from datetime import datetime
import requests


class PrintCallback(SwanKitCallback):
    """Basic callback for printing experiment status information."""

    def on_init(
        self, proj_name: str, workspace: str, logdir: Optional[str] = None, **kwargs
    ):
        """Called when experiment initialization completes."""
        print(f"🚀 My callback: on_init: {proj_name}, {workspace}, {logdir}, {kwargs}")

    def on_stop(self, error: Optional[str] = None):
        """Called when experiment stops or encounters error."""
        status = f"with error: {error}" if error else "successfully"
        print(f"🚀 My callback: Experiment stopped {status}")

    def __str__(self):
        return "PrintCallback"


class EmailCallback(SwanKitCallback):
    """Email notification callback with bilingual support."""

    DEFAULT_TEMPLATES = {
        "en": {
            "subject_success": "SwanLab | Your experiment completed successfully",
            "subject_error": "SwanLab | Your experiment encountered an error",
            "body_success": "Your SwanLab experiment has completed successfully.",
            "body_error": "Your SwanLab experiment encountered an error: {error}",
            "link_text": "\nExperiment Link: {link}",
        },
        "zh": {
            "subject_success": "SwanLab | 您的实验已成功完成",
            "subject_error": "SwanLab | 您的实验遇到错误",
            "body_success": "您的 SwanLab 实验已成功完成。",
            "body_error": "您的 SwanLab 实验遇到错误: {error}",
            "link_text": "\n实验链接: {link}",
        },
    }

    def __init__(
        self,
        sender_email: str,
        receiver_email: str,
        password: str,
        smtp_server: str = "smtp.gmail.com",
        port: int = 587,
        language: str = "en",
    ):
        """
        Initialize email callback configuration.

        :param sender_email: SMTP account email address
        :param receiver_email: Recipient email address
        :param password: SMTP account password
        :param smtp_server: SMTP server address
        :param port: SMTP server port
        :param language: Email content language (en/zh)
        """
        self.sender_email = sender_email
        self.receiver_email = receiver_email
        self.password = password
        self.smtp_server = smtp_server
        self.port = port
        self.language = language

    def _create_email_content(self, error: Optional[str] = None) -> Dict[str, str]:
        """Generate bilingual email content based on experiment status."""
        templates = self.DEFAULT_TEMPLATES[self.language]

        # Determine email subject and body based on error status
        if error:
            subject = templates["subject_error"]
            body = templates["body_error"].format(error=error)
        else:
            subject = templates["subject_success"]
            body = templates["body_success"]

        # Add experiment link if running in cloud mode
        exp_link = swanlab.get_url()
        if exp_link:
            body += templates["link_text"].format(link=exp_link)

        return {"subject": subject, "body": body}

    def _send_email(self, content: Dict[str, str]) -> None:
        """Handle SMTP connection and email sending."""
        try:
            # Create email message
            msg = MIMEMultipart()
            msg["From"] = self.sender_email
            msg["To"] = self.receiver_email
            msg["Subject"] = content["subject"]
            msg.attach(MIMEText(content["body"], "plain", "utf-8"))
            # Establish secure connection
            with smtplib.SMTP(self.smtp_server, self.port) as server:
                server.starttls()
                server.login(self.sender_email, self.password)
                server.send_message(msg)
                print("✅ Email sent successfully!")

        except smtplib.SMTPException as e:
            print(f"❌ Email sending failed: {str(e)}")

    def on_stop(self, error: Optional[str] = None) -> None:
        """Trigger email notification when experiment stops."""
        print("📧 Preparing email notification...")
        email_content = self._create_email_content(error)
        self._send_email(email_content)

    def __str__(self):
        return f"EmailCallback({self.receiver_email})"


class LarkBot:
    def __init__(self, webhook_url: str, secret: Optional[str] = None):
        self.webhook_url = webhook_url
        self.secret = secret

    def gen_sign(self, timesteamp: int) -> str:
        """
        docs: https://open.feishu.cn/document/client-docs/bot-v3/add-custom-bot?lang=zh-CN#9ff32e8e
        如果用户配置了签名校验功能，需要使用此方法生成签名
        """
        if not self.secret:
            raise ValueError("secret is required")
        string_to_sign: str = f"{timesteamp}\n{self.secret}"
        hmac_code = hmac.new(
            string_to_sign.encode("utf-8"), digestmod=hashlib.sha256
        ).digest()
        return base64.b64encode(hmac_code).decode("utf-8")


class LarkCallback(SwanKitCallback):
    DEFAULT_TEMPLATES = {
        "en": {
            "title": "SwanLab Message Notification\n",
            "msg_success": "SwanLab | Your experiment completed successfully\n",
            "msg_error": "Your SwanLab experiment encountered an error: {error}\n",
            "link_text": "Experiment Link: \n{link}",
        },
        "zh": {
            "title": "SwanLab 消息通知\n",
            "msg_success": "SwanLab | 您的实验已成功完成\n",
            "msg_error": "您的 SwanLab 实验遇到错误: {error}\n",
            "link_text": "实验链接: \n{link}",
        },
    }
    """Lark notification callback with bilingual support."""

    def __init__(
        self, webhook_url: str, secret: Optional[str] = None, language: str = "zh"
    ):
        self.lark_bot = LarkBot(webhook_url, secret)
        self.language = language

    def _create_lark_content(self, error: Optional[str] = None) -> str:
        templates = self.DEFAULT_TEMPLATES[self.language]
        content_text: str = f"{templates['title']}"
        if error:
            content_text += templates["msg_error"].format(error=error)
        else:
            content_text += templates["msg_success"]
        exp_link = swanlab.get_url()
        if exp_link:
            content_text += templates["link_text"].format(link=exp_link)
        return content_text

    def _send_lark_msg(self, content: str) -> None:
        timestamp: int = int(datetime.now().timestamp())
        json_params: Dict[str, Any] = {
            "timestamp": timestamp,
            "msg_type": "text",
            "content": {"text": content},
        }
        if self.lark_bot.secret:
            json_params["sign"] = self.lark_bot.gen_sign(timestamp)
        resp = requests.post(self.lark_bot.webhook_url, json=json_params)
        resp.raise_for_status()
        result: Dict[str, Any] = resp.json()
        if result.get("code") and result["code"] != 0:
            print(f"❌ LarkBot sending failed: {result.get('msg')}")
            return
        print("✅ LarkBot sending successfully")
        return

    def on_stop(self, error: Optional[str] = None) -> None:
        print("🤖 Preparing LarkBot notification...")
        lark_content = self._create_lark_content(error)
        self._send_lark_msg(lark_content)

    def __str__(self):
        return f"LarkBotCallback({self.lark_bot.webhook_url})"
