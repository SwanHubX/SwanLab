from datetime import datetime
import requests
from typing import Dict, Any, Optional
import hmac
import hashlib
import swanlab
import base64
from swankit.callback import SwanKitCallback


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


class LarkBotCallback(SwanKitCallback):
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
        self, webhook_url: str, secret: Optional[str] = None, language: str = "zh"
    ):
        super().__init__()
        self.language = language
        self.lark_bot = LarkBot(webhook_url, secret, language)

    def _create_lark_content(self, error: Optional[str] = None) -> Dict[str, Any]:
        templates = self.DEFAULT_TEMPLATES[self.language]
        content: Dict[str, Any] = {"post": {}}
        lang_key = "zh_cn" if self.language == "zh" else "en_us"
        content["post"][lang_key] = {}
        json_data = content["post"][lang_key]

        if error:
            json_data["title"] = templates["subject_error"]
            json_data["content"] = []
            json_data["content"].append(
                {"tag": "text", "text": templates["body_error"].format(error=error)}
            )
        else:
            json_data["title"] = templates["subject_success"]
            json_data["content"] = []
            json_data["content"].append(
                {"tag": "text", "text": templates["body_success"]}
            )
        exp_link = swanlab.get_url()
        if exp_link:
            json_data["content"].append(
                {
                    "tag": "a",
                    "text": templates["link_text"].format(link=exp_link),
                    "href": exp_link,
                }
            )

        return content

    def _send_lark_msg(self, content: Dict[str, str]) -> None:
        timestamp: int = int(datetime.now().timestamp())
        json_params: Dict[str, Any] = {
            "timestamp": timestamp,
            "msg_type": "post",
            "content": content,
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
