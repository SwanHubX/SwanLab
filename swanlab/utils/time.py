from datetime import datetime


def create_time():
    """获取当前时间"""
    return datetime.utcnow().isoformat() + "Z"

