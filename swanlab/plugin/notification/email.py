import html
import smtplib
import ssl
from email.message import EmailMessage
from typing import Optional

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe

_EMAIL_TEMPLATES = {
    "en": {
        "subject_success": "SwanLab | Experiment Completed Successfully",
        "subject_error": "SwanLab | Experiment Encountered an Error",
        "status_success": "✅ Experiment Completed Successfully",
        "status_error": "❌ Experiment Failed",
        "error_label": "Error Details",
        "project_label": "Project",
        "workspace_label": "Workspace",
        "exp_name_label": "Experiment",
        "desc_label": "Description",
        "view_btn": "View Experiment",
    },
    "zh": {
        "subject_success": "SwanLab | 实验已成功完成",
        "subject_error": "SwanLab | 实验遇到错误",
        "status_success": "✅ 实验已成功完成",
        "status_error": "❌ 实验失败",
        "error_label": "错误详情",
        "project_label": "项目",
        "workspace_label": "工作区",
        "exp_name_label": "实验名",
        "desc_label": "描述",
        "view_btn": "查看实验",
    },
}

_HTML_BODY = """\
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"></head>
<body style="margin:0;padding:32px;background:#f5f5f7;font-family:-apple-system,BlinkMacSystemFont,'SF Pro Text','SF Pro Display',system-ui,sans-serif;">
<div style="max-width:540px;margin:0 auto;background:#fff;border-radius:18px;overflow:hidden;">
  <div style="border-top:3px solid {status_accent};padding:32px 32px 8px;text-align:center;">
    <h1 style="color:#1d1d1f;margin:0;font-size:21px;font-weight:600;letter-spacing:0.231px;line-height:1.19;">{status_title}</h1>
  </div>
  <div style="padding:24px 32px 32px;">
    {error_section}
    <table style="width:100%;border-collapse:collapse;font-size:14px;line-height:1.43;">
      {metadata_rows}
    </table>
    {link_section}
  </div>
</div>
</body>
</html>"""


class EmailCallback(NotificationCallback):
    """Email notification callback with HTML formatting.

    Supports both STARTTLS (port 587) and direct SSL (port 465).

    Usage::

        from swanlab.plugin import EmailCallback

        swanlab.init(
            callbacks=[EmailCallback(
                sender_email="xxx@gmail.com",
                receiver_email="yyy@gmail.com",
                password="app-password",
                smtp_server="smtp.gmail.com",
                port=587,
            )]
        )
    """

    def __init__(
        self,
        sender_email: str,
        receiver_email: str,
        password: str,
        smtp_server: str = "smtp.gmail.com",
        port: int = 587,
        language: str = "zh",
    ):
        super().__init__(language=language)
        self._sender_email = sender_email
        self._receiver_email = receiver_email
        self._password = password
        self._smtp_server = smtp_server
        self._port = port

    def _build_subject(self, error: Optional[str]) -> str:
        tpl = _EMAIL_TEMPLATES.get(self.language, _EMAIL_TEMPLATES["en"])
        return tpl["subject_error"] if error else tpl["subject_success"]

    def _build_html(self, error: Optional[str]) -> str:
        tpl = _EMAIL_TEMPLATES.get(self.language, _EMAIL_TEMPLATES["en"])
        has_error = error is not None
        status_accent = "#ff3b30" if has_error else "#30d158"
        status_title = tpl["status_error"] if has_error else tpl["status_success"]

        error_section = ""
        if has_error:
            escaped_error = html.escape(error)
            error_section = (
                '<div style="padding:12px 0 16px;border-bottom:1px solid #f0f0f0;'
                f'margin-bottom:20px;font-size:14px;color:#1d1d1f;line-height:1.47;">'
                f'<span style="font-weight:600;color:#ff3b30;">{tpl["error_label"]}</span><br>'
                f"{escaped_error}</div>"
            )

        rows = ""
        if self._settings:
            meta = [
                (tpl["project_label"], self._settings.project.name or "-"),
                (tpl["workspace_label"], self._settings.project.workspace or "-"),
                (tpl["exp_name_label"], self._settings.experiment.name or "-"),
                (tpl["desc_label"], self._settings.experiment.description or "-"),
            ]
            for i, (label, value) in enumerate(meta):
                border = "border-bottom:1px solid #f0f0f0;" if i < len(meta) - 1 else ""
                rows += (
                    f'<tr style="{border}">'
                    f'<td style="padding:10px 0;color:#7a7a7a;width:90px;vertical-align:top;">'
                    f"{label}</td>"
                    f'<td style="padding:10px 0;color:#1d1d1f;">'
                    f"{html.escape(value)}</td></tr>"
                )

        link_section = ""
        url = self._get_url()
        if url:
            link_section = (
                '<div style="margin-top:28px;text-align:center;">'
                f'<a href="{url}" style="display:inline-block;padding:11px 22px;background:#0066cc;'
                f"color:#fff;text-decoration:none;border-radius:9999px;font-size:14px;"
                f'letter-spacing:-0.224px;">'
                f"{tpl['view_btn']}</a></div>"
            )

        return _HTML_BODY.format(
            status_accent=status_accent,
            status_title=status_title,
            error_section=error_section,
            metadata_rows=rows,
            link_section=link_section,
        )

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        subject = self._build_subject(error)

        msg = EmailMessage()
        msg["From"] = self._sender_email
        msg["To"] = self._receiver_email
        msg["Subject"] = subject
        msg.set_content(self._build_content(state, error))
        msg.add_alternative(self._build_html(error), subtype="html")

        with safe.block(smtplib.SMTPException, message="❌ EmailBot sending failed"):
            ctx = ssl.create_default_context()
            if self._port == 465:
                with smtplib.SMTP_SSL(self._smtp_server, self._port, context=ctx, timeout=15) as server:
                    server.login(self._sender_email, self._password)
                    server.send_message(msg)
            else:
                with smtplib.SMTP(self._smtp_server, self._port, timeout=15) as server:
                    server.starttls(context=ctx)
                    server.login(self._sender_email, self._password)
                    server.send_message(msg)
            console.info("✅ EmailBot notification sent successfully")
