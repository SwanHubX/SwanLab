from datetime import datetime
import pytz


def create_time():
    """获取当前时间"""
    # 获取当前时间
    current_time = datetime.utcnow()
    # 获取UTC时区对象
    utc_timezone = pytz.utc
    # 将当前时间转换为带有时区信息的时间
    current_time_with_timezone = utc_timezone.localize(current_time)
    # 返回带有时区信息的当前时间
    return current_time_with_timezone.isoformat()
