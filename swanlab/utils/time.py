from datetime import datetime
import pytz


def create_time():
    """获取当前时间"""
    return datetime.utcnow().isoformat() + "Z"

