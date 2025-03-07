from swankit.callback import SwanKitCallback
import swanlab


class PrintCallback(SwanKitCallback):
    def on_init(self, proj_name: str, workspace: str, logdir: str = None, **kwargs):
        print(f"🚀my callback: on_init: {proj_name}, {workspace}, {logdir}, {kwargs}")
      
    def on_stop(self, error: str = None):
        print(f"🚀my callback: on_stop: {error}")
      
    def __str__(self):
        return "PrintCallback"


class EmailCallback(SwanKitCallback):
    def __init__(self, sender_email: str, receiver_email: str, password: str, smtp_server: str = "smtp.gmail.com", port: int = 587):
        self.sender_email = sender_email
        self.receiver_email = receiver_email
        self.password = password
        self.smtp_server = smtp_server
        self.port = port
    
    def _create_template(self, error: str = None, mode: str = "cloud", exp_link: str = None):
        # 根据实验状态和模式生成邮件模板
        status = "encountered an error" if error else "completed successfully"
        
        subject = f"SwanLab | Your experiment {status}"
        body = f"Your SwanLab experiment {status}."
        
        # 仅在云模式下添加实验链接
        if mode == "cloud":
            body += f"\nExperiment Link: {exp_link}"
        
        return subject, body
      
    def on_stop(self, error: str = None):
        print(f"📧EmailCallback Launch")
        import smtplib
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart
        
        # 创建邮件内容
        message = MIMEMultipart()
        message["From"] = self.sender_email
        message["To"] = self.receiver_email
        
        print(swanlab.get_url())
        
        if swanlab.get_url() is None:
            mode = "local"
            exp_link = None
        else:
            mode = "cloud"
            exp_link = swanlab.get_url()
        
        subject, body = self._create_template(error, mode, exp_link)
            
        try:
            # 创建邮件对象
            msg = MIMEMultipart()
            msg['From'] = self.sender_email
            msg['To'] = self.receiver_email
            msg['Subject'] = subject
            msg.attach(MIMEText(body, 'plain', 'utf-8'))

            # 连接SMTP服务器并发送邮件
            with smtplib.SMTP_SSL(self.smtp_server, self.port) as server:
                print(f"🔑Login: {self.sender_email}, {self.password}")
                server.login(self.sender_email, self.password)
                server.sendmail(self.sender_email, self.receiver_email, msg.as_string())
                print("Email sent successfully!")

        except Exception as e:
            print(f"Email sending failed: {e}")
    
    def __str__(self):
        return "EmailCallback"