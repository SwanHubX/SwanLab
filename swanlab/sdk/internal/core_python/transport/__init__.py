"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/16
@description: 事件驱动的 Record 上传模块

数据流：put() → RecordBuffer 攒批去重 → Dispatch 分组分片 → Sender 上传 → UploadResult。

上传结果统计：
- confirmed_records：本轮已终结，不再重试（上传成功，或确定性失败）；
- remaining_records：交回 worker 退避重试；
- failures：失败诊断（message / status / code / exception），不参与重试决策。

重试边界由 sender.utils.is_retryable_error 统一判定：
- 网络异常与 5xx 保留重试；
- 本地数据不可用、协议错误、save 文件上传失败等无法自愈的失败终结记录。
- 后端服务长时间不可用时，record 留在进程内缓冲，不阻塞训练线程。

进度统计：
- total_records 为入队去重后的记录数，
- uploaded_records 累计 confirmed_records，代表"已确认"而非"上传成功"，差值即待重试数；
- finish 结束信号只判断 pending 是否排空，不由统计字段推断。
"""

from .worker import Transport

__all__ = ["Transport"]
