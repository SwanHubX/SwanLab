"""
Notification plugin for SwanLab.
Used for sending notifications to users.
"""

import base64
import hashlib
import hmac
import smtplib
from abc import ABC, abstractmethod
from datetime import datetime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from typing import Optional, Dict, Any, Tuple, Literal

import requests
import json
import swanlab
from swanlab.toolkit import SwanKitCallback


class PrintCallback(SwanKitCallback):
    """Basic callback for printing experiment status information."""

    def on_init(self, proj_name: str, workspace: str, logdir: Optional[str] = None, *args, **kwargs):
        """Called when experiment initialization completes."""
        print(f"🚀 My callback: on_init: {proj_name}, {workspace}, {logdir}, {kwargs}")

    def on_stop(self, error: Optional[str] = None, *args, **kwargs):
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

    def _create_email_content(self, error: Optional[str] = None) -> Tuple[str, str]:
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

    def on_init(self, proj_name: str, workspace: str, logdir: Optional[str] = None, *args, **kwargs):
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_stop(self, error: Optional[str] = None, *args, **kwargs):
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

    def on_init(self, proj_name: str, workspace: str, logdir: Optional[str] = None, *args, **kwargs):
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_stop(self, error: Optional[str] = None, *args, **kwargs):
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
        if resp.status_code not in [200, 204]:
            print(f"❌ DiscordBot sending failed: {resp.text}")
            return
        print("✅ DiscordBot sending successfully")

    def __str__(self):
        return "DiscordBotCallback"


class SlackBot:
    """
    Slack notification callback with bilingual support.
    docs: https://api.slack.com/messaging/webhooks
    """

    def __init__(self, webhook_url: str):
        self.webhook_url = webhook_url


class SlackCallback(WebhookCallback):
    """Slack notification callback with bilingual support."""

    def _init_bot(self) -> None:
        self.bot = SlackBot(self.webhook_url)

    def send_msg(self, content: str) -> None:
        data = {
            "text": content,
        }
        resp = requests.post(self.bot.webhook_url, json=data)
        resp.raise_for_status()
        if resp.status_code not in [200, 204]:
            print(f"❌ SlackBot sending failed: {resp.text}")
            return
        print("✅ SlackBot sending successfully")

    def __str__(self):
        return "SlackBotCallback"


class BarkCallback(SwanKitCallback):
    """使用Bark对IOS设备进行通知"""

    def __str__(self) -> str:
        return "BarkCallback"

    DEFAULT_TEMPLATES = {
        'en': {
            "subtitle_success": "Your experiment completed successfully",
            "subtitle_error": "Your experiment encountered an error: {error}",
            "link_text": "Project: {project}\nWorkspace: {workspace}\nName: {exp_name}\nDescription: {description}",
            "offline_text": "Offline use of the project."
        },
        "zh": {
            "subtitle_success": "您的实验已成功完成",
            "subtitle_error": "您的实验遇到错误: {error}",
            "link_text": "项目: {project}\n工作区: {workspace}\n实验名: {exp_name}\n描述: {description}",
            "offline_text": "项目离线使用。"
        },
    }
    def __init__(
            self,
            key: str,
            url: str = 'https://api.day.app',
            title: str = 'SwanLab',
            bark_level: Literal['critical', 'active', 'timeSensitive', 'passive'] = 'active',
            icon: str = 'https://www.swanlab.cn/icon.png',
            group: Literal['exp_name', 'project', 'workspace', None] = None,
            click_jump: bool = True,
            language: str = 'zh',
    ):
        """
        初始化 Bark callback 配置
        详细内容见 https://bark.day.app/#/tutorial?id=%e8%af%b7%e6%b1%82%e5%8f%82%e6%95%b0
        :param key: bark中的device key
        :param url: bark推送链接
        :param title: 推送的通知标题，默认为SwanLab
        :param bark_level: bark推送等级
        :param icon: bark推送图标，默认为SwanLab图标
        :param group: bark推送分组，默认为None为不分组
        :param click_jump: 是否点击链接跳转，默认为True
        :param language: 推送语言
        """
        self.key = key
        self.url = url
        self.title = title
        if bark_level not in ['critical', 'active', 'timeSensitive', 'passive']:
            raise ValueError(f'Invalid bark_level {bark_level}')
        self.bark_level = bark_level
        self.icon = icon
        self.group = group
        self.click_jump = click_jump
        self.language = language


    def _create_notification_message(self, error: Optional[str] = None) -> Dict[str, Optional[str]]:
        """根据实验状态，生成对应的通知消息"""
        templates = self.DEFAULT_TEMPLATES[self.language]
        if error:
            subtitle = templates['subtitle_error'].format(error=error)
        else:
            subtitle = templates['subtitle_success']

        exp_link = swanlab.get_url()
        if exp_link:
            body = templates['link_text'].format(
                project=self.project,
                workspace=self.workspace,
                exp_name=self.exp_name,
                description=self.description
            )
            if self.group == 'exp_name':
                group = self.exp_name
            elif self.group == 'project':
                group = self.project
            elif self.group == 'workspace':
                group = self.workspace
            else:
                group = None
        else:
            body = templates['offline_text']
            group = None

        # 构建post字典
        data = {
            'title': self.title,
            'subtitle': subtitle,
            'body': body,
            'level': self.bark_level,
            'icon': self.icon,
            'group': group,
            'device_key': self.key,
            'url': exp_link if self.click_jump else None,
        }
        data = {k: v for k, v in data.items() if v is not None}
        return data


    def send_notification(self, data: dict):
        """发送通知，也可以直接构建数据字典进行消息发送"""
        if self.url.endswith('/'):
            url = self.url + 'push'
        else:
            url = self.url + '/push'

        headers = {
            'Content-Type': 'application/json; charset=utf-8',
        }

        resp = requests.post(
            url=url,
            headers=headers,
            data=json.dumps(data, ensure_ascii=False).encode('utf-8'),
        )
        resp.raise_for_status()
        result: Dict[str, Any] = resp.json()
        if result.get("errcode") and result["errcode"] != 0:
            print(f"❌ Bark sending failed: {result.get('errmsg')}")
            return
        print("✅ Bark sending successfully")


    def on_init(self, proj_name: str, workspace: str, public: Optional[bool] = None, logdir: Optional[str] = None, *args, **kwargs):
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_stop(self, error: Optional[str] = None, *args, **kwargs):
        content = self._create_notification_message(error)
        self.send_notification(content)


class TelegramBot:
    """
    Telegram Bot notification helper.
    docs: https://core.telegram.org/bots/api#sendmessage
    """

    def __init__(self, bot_token: str, chat_id: str):
        """
        Initialize Telegram Bot.

        :param bot_token: Telegram Bot API token (get from @BotFather)
        :param chat_id: Target chat ID (can be user ID, group ID, or channel username)
        """
        self.bot_token = bot_token
        self.chat_id = chat_id
        self.api_base = f"https://api.telegram.org/bot{bot_token}"

    def send_message(self, text: str, parse_mode: str = "HTML") -> Dict[str, Any]:
        """
        Send a message via Telegram Bot API.

        :param text: Message text
        :param parse_mode: Parse mode for formatting (HTML or Markdown)
        :return: API response
        """
        url = f"{self.api_base}/sendMessage"
        data = {
            "chat_id": self.chat_id,
            "text": text,
            "parse_mode": parse_mode,
        }
        resp = requests.post(url, json=data)
        resp.raise_for_status()
        return resp.json()


class TelegramCallback(SwanKitCallback):
    """
    Telegram notification callback with bilingual support.
    Send notifications to Telegram chat when experiment starts/stops.

    Usage:
        1. Create a bot via @BotFather and get the bot token
        2. Get your chat ID (send /start to @userinfobot)
        3. Initialize the callback:

        ```python
        from swanlab.plugin import TelegramCallback

        telegram = TelegramCallback(
            bot_token="YOUR_BOT_TOKEN",
            chat_id="YOUR_CHAT_ID",
            language="zh"  # or "en"
        )
        swanlab.init(callbacks=[telegram])
        ```
    """

    DEFAULT_TEMPLATES = {
        "en": {
            "title": "🧪 <b>SwanLab Notification</b>\n\n",
            "msg_start": "🚀 Experiment started\n",
            "msg_success": "✅ Experiment completed successfully\n",
            "msg_error": "❌ Experiment failed: {error}\n",
            "link_text": (
                "<b>Project:</b> {project}\n"
                "<b>Workspace:</b> {workspace}\n"
                "<b>Name:</b> {exp_name}\n"
                "<b>Description:</b> {description}\n"
                "<b>Link:</b> {link}"
            ),
            "offline_text": "📴 Running in offline mode",
        },
        "zh": {
            "title": "🧪 <b>SwanLab 消息通知</b>\n\n",
            "msg_start": "🚀 实验已开始\n",
            "msg_success": "✅ 实验已成功完成\n",
            "msg_error": "❌ 实验遇到错误: {error}\n",
            "link_text": (
                "<b>项目:</b> {project}\n"
                "<b>工作区:</b> {workspace}\n"
                "<b>实验名:</b> {exp_name}\n"
                "<b>描述:</b> {description}\n"
                "<b>链接:</b> {link}"
            ),
            "offline_text": "📴 离线模式运行中",
        },
    }

    def __init__(
        self,
        bot_token: str,
        chat_id: str,
        language: str = "zh",
        notify_on_start: bool = False,
    ):
        """
        Initialize Telegram callback configuration.

        :param bot_token: Telegram Bot API token (get from @BotFather)
        :param chat_id: Target chat ID (user ID, group ID, or @channel_username)
        :param language: Notification language (en/zh)
        :param notify_on_start: Whether to send notification when experiment starts
        """
        self.bot = TelegramBot(bot_token, chat_id)
        self.language = language
        self.notify_on_start = notify_on_start

    def _create_content(self, event: str = "stop", error: Optional[str] = None) -> str:
        """
        Create notification content based on event type.

        :param event: Event type ("start" or "stop")
        :param error: Error message if experiment failed
        :return: Formatted message text
        """
        templates = self.DEFAULT_TEMPLATES[self.language]
        content = templates["title"]

        if event == "start":
            content += templates["msg_start"]
        elif error:
            content += templates["msg_error"].format(error=error)
        else:
            content += templates["msg_success"]

        exp_link = swanlab.get_url()
        if exp_link:
            content += templates["link_text"].format(
                project=self.project,
                workspace=self.workspace,
                exp_name=self.exp_name,
                description=self.description or "N/A",
                link=exp_link,
            )
        else:
            content += templates["offline_text"]

        return content

    def send_msg(self, content: str) -> None:
        """Send message via Telegram Bot."""
        try:
            result = self.bot.send_message(content)
            if result.get("ok"):
                print("✅ Telegram message sent successfully")
            else:
                print(f"❌ Telegram sending failed: {result.get('description')}")
        except requests.RequestException as e:
            print(f"❌ Telegram sending failed: {str(e)}")

    def on_init(self, proj_name: str, workspace: str, public: Optional[bool] = None, logdir: Optional[str] = None, *args, **kwargs):
        """Called when experiment initialization completes."""
        self.project = proj_name
        self.workspace = workspace

    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        """Called before experiment initialization."""
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description

    def on_run(self, *args, **kwargs):
        """Called when experiment starts running."""
        if self.notify_on_start:
            print("📱 Preparing Telegram start notification...")
            content = self._create_content(event="start")
            self.send_msg(content)

    def on_stop(self, error: Optional[str] = None, *args, **kwargs):
        """Called when experiment stops."""
        print("📱 Preparing Telegram notification...")
        content = self._create_content(event="stop", error=error)
        self.send_msg(content)

    def __str__(self) -> str:
        return "TelegramCallback"
