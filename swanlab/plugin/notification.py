"""
Notification plugin for SwanLab.
Used for sending notifications to users.
"""

from swankit.callback import SwanKitCallback
import swanlab
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Optional, Dict, Any, Tuple
import hmac
import hashlib
import base64
from datetime import datetime
import requests
from abc import ABC, abstractmethod


class PrintCallback(SwanKitCallback):
    """Basic callback for printing experiment status information."""

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs):
        """Called when experiment initialization completes."""
        print(f"🚀 My callback: on_init: {proj_name}, {workspace}, {logdir}, {kwargs}")

    def on_stop(self, error: str = None, *args, **kwargs):
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
            "body_success": "Your SwanLab experiment has completed successfully.\n",
            "body_error": "Your SwanLab experiment encountered an error: {error}\n",
            "link_text": "Project: {project}\nWorkspace: {workspace}\nName: {exp_name}\nDescription: {description}\nExperiment Link: {link}",
        },
        "zh": {
            "subject_success": "SwanLab | 您的实验已成功完成",
            "subject_error": "SwanLab | 您的实验遇到错误",
            "body_success": "您的 SwanLab 实验已成功完成。\n",
            "body_error": "您的 SwanLab 实验遇到错误: {error}\n",
            "link_text": "项目: {project}\n工作区: {workspace}\n实验名: {exp_name}\n描述: {description}\n实验链接: {link}",
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
            body += templates["link_text"].format(
                project=self.project,
                workspace=self.workspace,
                exp_name=self.exp_name,
                description=self.description,
                link=exp_link,
            )

        return subject, body

    def send_email(self, subject: str, body: str) -> None:
        """Handle SMTP connection and email sending."""
        try:
            # Create email message
            msg = MIMEMultipart()
            msg["From"] = self.sender_email
            msg["To"] = self.receiver_email
            msg["Subject"] = subject
            msg.attach(MIMEText(body, "plain", "utf-8"))
            # Establish secure connection
            with smtplib.SMTP(self.smtp_server, self.port) as server:
                server.starttls()
                server.login(self.sender_email, self.password)
                server.send_message(msg)
                print("✅ Email sent successfully!")

        except smtplib.SMTPException as e:
            print(f"❌ Email sending failed: {str(e)}")

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs):
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_stop(self, error: str = None, *args, **kwargs):
        """Trigger email notification when experiment stops."""
        print("📧 Preparing email notification...")
        subject, body = self._create_email_content(error)
        self.send_email(subject, body)

    def __str__(self):
        return "EmailCallback"


WEBHOOK_DEFAULT_TEMPLATES = {
    "en": {
        "title": "SwanLab Message Notification\n",
        "msg_success": "SwanLab | Your experiment completed successfully\n",
        "msg_error": "Your SwanLab experiment encountered an error: {error}\n",
        "link_text": "Project: {project}\nWorkspace: {workspace}\nName: {exp_name}\nDescription: {description}\nExperiment Link: {link}",
    },
    "zh": {
        "title": "SwanLab 消息通知\n",
        "msg_success": "SwanLab | 您的实验已成功完成\n",
        "msg_error": "您的 SwanLab 实验遇到错误: {error}\n",
        "link_text": "项目: {project}\n工作区: {workspace}\n实验名: {exp_name}\n描述: {description}\n实验链接: {link}",
    },
}


class WebhookCallback(SwanKitCallback, ABC):
    """抽象基类，用于各种webhook通知回调"""

    def __init__(self, webhook_url: str, secret: Optional[str] = None, language: str = "zh"):
        self.webhook_url = webhook_url
        self.secret = secret
        self.language = language
        self._init_bot()

    @abstractmethod
    def _init_bot(self) -> None:
        """初始化机器人实例"""
        pass

    def _create_content(self, error: Optional[str] = None) -> str:
        """创建通知内容"""
        templates = WEBHOOK_DEFAULT_TEMPLATES[self.language]
        content_text: str = f"{templates['title']}"
        if error:
            content_text += templates["msg_error"].format(error=error)
        else:
            content_text += templates["msg_success"]
        exp_link = swanlab.get_url()
        if exp_link:
            content_text += templates["link_text"].format(
                project=self.project,
                workspace=self.workspace,
                exp_name=self.exp_name,
                description=self.description,
                link=exp_link,
            )
        return content_text

    @abstractmethod
    def send_msg(self, content: str) -> None:
        """发送消息的具体实现"""
        pass

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs):
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_stop(self, error: str = None, *args, **kwargs):
        print(f"🤖 Preparing {self.__class__.__name__} notification...")
        content = self._create_content(error)
        self.send_msg(content)


class LarkBot:
    def __init__(self, webhook_url: str, secret: Optional[str] = None):
        self.webhook_url = webhook_url
        self.secret = secret

    def gen_sign(self, timesteamp: int) -> str:
        """
        docs: https://open.feishu.cn/document/client-docs/bot-v3/add-custom-bot?lang=zh-CN#9ff32e8e
        If the user has configured the signature verification function, this method is required to generate the signature
        """
        if not self.secret:
            raise ValueError("secret is required")
        string_to_sign: str = f"{timesteamp}\n{self.secret}"
        hmac_code = hmac.new(string_to_sign.encode("utf-8"), digestmod=hashlib.sha256).digest()
        return base64.b64encode(hmac_code).decode("utf-8")


class LarkCallback(WebhookCallback):
    """Lark notification callback with bilingual support."""

    def _init_bot(self) -> None:
        self.bot = LarkBot(self.webhook_url, self.secret)

    def send_msg(self, content: str) -> None:
        timestamp: int = int(datetime.now().timestamp())
        json_params: Dict[str, Any] = {
            "timestamp": timestamp,
            "msg_type": "text",
            "content": {"text": content},
        }
        if self.bot.secret:
            json_params["sign"] = self.bot.gen_sign(timestamp)
        resp = requests.post(self.bot.webhook_url, json=json_params)
        resp.raise_for_status()
        result: Dict[str, Any] = resp.json()
        if result.get("code") and result["code"] != 0:
            print(f"❌ LarkBot sending failed: {result.get('msg')}")
            return
        print("✅ LarkBot sending successfully")

    def __str__(self):
        return "LarkBotCallback"


class DingTalkBot:
    """DingTalk notification callback with bilingual support."""

    def __init__(self, webhook_url: str, secret: Optional[str] = None):
        self.webhook_url = webhook_url
        self.secret = secret
        self.start_time = datetime.now().timestamp()  # 加签时，请求时间戳与请求时间不能超过1小时，用于定时更新签名
        if self.secret is not None and self.secret.startswith("SEC"):
            self.update_webhook_url()

    def _check_sign(self) -> None:
        """
        DingTalk 中的签名和时间戳需要放在 url 中，且签名中的时间戳与请求时间戳的差值不能超过 1 小时
        """
        ts_now = datetime.now().timestamp()
        if ts_now - self.start_time >= 3600 and self.secret is not None and self.secret.startswith("SEC"):
            self.start_time = ts_now
            self.update_webhook_url()

    def update_webhook_url(self) -> None:
        """
        DingTalk 中的 webhook url 需要更新签名和时间戳
        """
        """
        docs: https://open.dingtalk.com/document/
        DingTalk 中的签名和时间戳需要放在 url 中，且签名中的时间戳与请求时间戳的差值不能超过 1 小时
        """
        if (not self.secret) or (not self.secret.startswith("SEC")):
            raise ValueError("secret is invalid")
        timestamp = round(self.start_time * 1000)
        # 根据钉钉文档，正确的签名格式应为 timestamp\nsecret
        string_to_sign = f"{timestamp}\n{self.secret}"
        hmac_code = hmac.new(self.secret.encode(), string_to_sign.encode(), digestmod=hashlib.sha256).digest()
        sign = base64.b64encode(hmac_code).decode("utf-8")

        # 提取基础 URL（移除旧的时间戳和签名参数）
        base_url = self.webhook_url.split("&timestamp=")[0] if "&timestamp=" in self.webhook_url else self.webhook_url

        # 添加分隔符（? 或 &）
        separator = "?" if "?" not in base_url else "&"
        if separator == "&" and base_url.endswith("&"):
            separator = ""

        # 构建新的 URL
        self.webhook_url = f"{base_url}{separator}timestamp={timestamp}&sign={sign}"


class DingTalkCallback(WebhookCallback):
    """DingTalk notification callback with bilingual support."""

    def _init_bot(self) -> None:
        self.bot = DingTalkBot(self.webhook_url, self.secret)

    def send_msg(self, content: str) -> None:
        self.bot._check_sign()

        data = {
            "msgtype": "text",
            "text": {"content": content},
        }
        resp = requests.post(self.bot.webhook_url, json=data)
        resp.raise_for_status()
        result: Dict[str, Any] = resp.json()
        if result.get("errcode") and result["errcode"] != 0:
            print(f"❌ DingTalkBot sending failed: {result.get('errmsg')}")
            return
        print("✅ DingTalkBot sending successfully")

    def __str__(self):
        return "DingTalkBotCallback"


class WXWorkBot:
    """
    WxWork notification callback with bilingual support.
    docs: https://developer.work.weixin.qq.com/document/path/91770
    """

    def __init__(self, webhook_url: str):
        # WxWork does not need secret
        self.webhook_url = webhook_url


class WXWorkCallback(WebhookCallback):
    """WxWork notification callback with bilingual support."""

    def _init_bot(self) -> None:
        self.bot = WXWorkBot(self.webhook_url)

    def send_msg(self, content: str) -> None:
        data = {
            "msgtype": "text",
            "text": {"content": content},
        }
        resp = requests.post(self.bot.webhook_url, json=data)
        resp.raise_for_status()
        result: Dict[str, Any] = resp.json()
        if result.get("errcode") and result["errcode"] != 0:
            print(f"❌ WXWorkBot sending failed: {result.get('errmsg')}")
            return
        print("✅ WXWorkBot sending successfully")

    def __str__(self):
        return "WXWorkBotCallback"


class DiscordBot:
    """
    Discord notification callback with bilingual support.
    docs: https://discord.com/developers/docs/resources/webhook
    """

    def __init__(self, webhook_url: str):
        self.webhook_url = webhook_url


class DiscordCallback(WebhookCallback):
    """Discord notification callback with bilingual support."""

    def _init_bot(self) -> None:
        self.bot = DiscordBot(self.webhook_url)

    def send_msg(self, content: str) -> None:
        data = {
            "content": content,
        }
        resp = requests.post(self.bot.webhook_url, json=data)
        resp.raise_for_status()
        print("✅ DiscordBot sending successfully")

    def __str__(self):
        return "DiscordBotCallback"
