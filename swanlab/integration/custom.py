from swankit.callback import SwanKitCallback
import swanlab
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Optional, Dict

class PrintCallback(SwanKitCallback):
    """Basic callback for printing experiment status information."""
    
    def on_init(self, proj_name: str, workspace: str, logdir: Optional[str] = None, **kwargs):
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
            "link_text": "\nExperiment Link: {link}"
        },
        "zh": {
            "subject_success": "SwanLab | 您的实验已成功完成",
            "subject_error": "SwanLab | 您的实验遇到错误",
            "body_success": "您的 SwanLab 实验已成功完成。",
            "body_error": "您的 SwanLab 实验遇到错误: {error}",
            "link_text": "\n实验链接: {link}"
        }
    }

    def __init__(
        self,
        sender_email: str,
        receiver_email: str,
        password: str,
        smtp_server: str = "smtp.gmail.com",
        port: int = 587,
        language: str = "en"
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
            msg['From'] = self.sender_email
            msg['To'] = self.receiver_email
            msg['Subject'] = content["subject"]
            msg.attach(MIMEText(content["body"], 'plain', 'utf-8'))
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