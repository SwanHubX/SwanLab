"""
@author: cunyue
@file: score.py
@description: 评测分数的事件时间分配器

分数的版本语义由**显式毫秒事件时间**承载：House 以 (projectId, evalId, dataset, key)
定位一个槽位，按毫秒时间戳取最新版本（ReplacingMergeTree + argMax）。因此：

1. 同一个槽位的连续更新必须严格递增，包括同一毫秒内连续写入与时钟回拨；
2. 时间戳在记录创建时分配一次，之后随记录一起持久化与重放，绝不重新生成——
   House 的 200 只代表"可能已落库"，换时间戳重发会让重放版本覆盖合法旧版本；
3. 跨进程重启后仍要严格递增，否则重启后写入的一条会与重启前的版本撞在同一毫秒上，
   胜负未定义。这里从本地数据存储里回放已写入的分数记录来恢复水位。
"""

import json
import threading
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Tuple

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScoreRecord
from swanlab.sdk.internal.pkg import console

__all__ = ["ScoreClock", "score_slot"]

ScoreSlot = Tuple[str, str, str]

# House 侧毫秒是版本精度，亚毫秒没有定义赢家，因此水位以毫秒为唯一单位
_MILLIS_PER_SECOND = 1000
_MILLIS_PER_NANO = 1_000_000


def score_slot(score: ScoreRecord) -> ScoreSlot:
    """分数的逻辑槽位：同一槽位内的写入共享一个单调时间序列。"""
    return score.eval_id, score.dataset, score.key


class ScoreClock:
    """按槽位分配严格递增的毫秒事件时间。

    时间戳取值 = max(当前毫秒, 该槽位上次分配值 + 1ms)，因此：
    - 正常写入时贴近墙钟；
    - 同一毫秒内的多次更新、以及系统时钟回拨时，仍然严格递增。
    """

    def __init__(
        self,
        state_path: Optional[Path] = None,
        recover: Optional[Callable[[], Iterable[ScoreRecord]]] = None,
    ) -> None:
        # 水位落盘：重启后仍需对同一槽位严格递增，否则重启后的写入可能和重启前的
        # 版本撞在同一毫秒上（House 的毫秒是版本精度，同毫秒胜负未定义）。
        #
        # 状态文件只是快路径。文件丢失、损坏或写盘失败时，用 recover() 从已落盘的分数
        # 记录重建水位——记录本身是持久化的事实来源，账本只是它的缓存。
        self._state_path = state_path
        self._recover = recover
        self._watermarks: Dict[ScoreSlot, int] = {}
        self._lock = threading.Lock()
        self._degraded = False
        # 账本只可能落后于记录（写失败、被回滚、进程崩溃前未落盘），因此真正的水位
        # 必须在第一次分配前与已落盘记录取最大值，不能只信账本。
        self._merged = recover is None
        # 恢复失败时无法确定安全水位：宁可拒写该进程的所有分数，也不冒"复用旧版本号"的风险
        self._unsafe_reason: Optional[str] = None
        self._load()

    @staticmethod
    def _now_ms() -> int:
        ts = Timestamp()
        ts.GetCurrentTime()
        return ts.seconds * _MILLIS_PER_SECOND + ts.nanos // _MILLIS_PER_NANO

    def _load(self) -> None:
        if self._state_path is None:
            self._restore_from_store()
            return
        if not self._state_path.exists():
            self._restore_from_store()
            return
        try:
            raw = json.loads(self._state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            console.warning(
                f"Score version ledger {self._state_path} is unreadable; rebuilding it from persisted score records."
            )
            self._restore_from_store()
            return
        for slot_str, watermark in raw.items():
            parts = slot_str.split("\x00")
            if len(parts) != 3 or not isinstance(watermark, int):
                continue
            self._watermarks[(parts[0], parts[1], parts[2])] = watermark

    def _ensure_merged(self) -> None:
        """把已落盘记录的最大水位并入账本水位（每个进程只做一次）。

        账本可能落后：水位写盘失败时记录已落盘、账本没更新；重启后若只读账本，
        同一槽位会重新分配一个已经用过的毫秒，两个不同分值撞在同一版本上。
        """
        if self._merged:
            return
        self._merged = True
        if self._recover is None:
            return
        self._restore_from_store()

    def _restore_from_store(self) -> None:
        """从已落盘的分数记录重建水位（账本丢失/损坏时的兜底）。"""
        if self._recover is None:
            return
        try:
            records = list(self._recover())
        except Exception as e:  # noqa: BLE001 - 恢复失败必须可见，且要阻止不安全写入
            self._unsafe_reason = f"cannot recover score version watermarks from persisted records: {e}"
            console.error(self._unsafe_reason)
            return
        self.restore(records)

    def _save(self) -> None:
        if self._state_path is None:
            return
        payload = {"\x00".join(slot): watermark for slot, watermark in self._watermarks.items()}
        tmp = self._state_path.with_suffix(".tmp")
        try:
            tmp.write_text(json.dumps(payload), encoding="utf-8")
            tmp.replace(self._state_path)
            self._degraded = False
        except OSError as e:
            # 账本落盘失败不阻断训练：记录本身已经（将）写入数据存储，重启时由
            # _restore_from_store 重建水位。但要留痕，否则"静默降级"看起来像正常工作。
            if not self._degraded:
                console.warning(f"Failed to persist score version ledger to {self._state_path}: {e}")
            self._degraded = True

    def restore(self, records: Iterable[ScoreRecord]) -> None:
        """用已落盘的历史分数恢复各槽位水位（进程启动时调用）。"""
        with self._lock:
            for record in records:
                if not record.HasField("timestamp"):
                    continue
                slot = score_slot(record)
                assigned = record.timestamp.seconds * _MILLIS_PER_SECOND + (record.timestamp.nanos // _MILLIS_PER_NANO)
                if assigned > self._watermarks.get(slot, -1):
                    self._watermarks[slot] = assigned

    def next_timestamp(self, slot: ScoreSlot) -> Timestamp:
        """为槽位分配新的毫秒时间戳，保证严格大于该槽位此前任何一条。"""
        self._ensure_merged()
        if self._unsafe_reason is not None:
            # 不知道安全水位就不写：写一个可能重复的版本会静默覆盖别人的分值
            raise RuntimeError(f"refusing to assign a score version: {self._unsafe_reason}")
        with self._lock:
            assigned = max(self._now_ms(), self._watermarks.get(slot, -1) + 1)
            self._watermarks[slot] = assigned
        self._save()
        ts = Timestamp()
        ts.FromNanoseconds(assigned * _MILLIS_PER_NANO)
        return ts

    def watermark(self, slot: ScoreSlot) -> Optional[int]:
        """当前槽位水位（毫秒），未写过时为 None。测试与诊断用。"""
        with self._lock:
            return self._watermarks.get(slot)
