"""W&B 事务日志（``.wandb`` 文件）只读扫描器。

文件为带 W&B 文件头的 leveldb journal，分帧符合 LevelDB 格式。
此处声明 W&B 的文件头签名（``:W&B`` / 0xBEE1 / version 0），并对损坏与截断加严：
转换场景下静默丢弃数据不可接受。
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path
from typing import IO, Optional, Tuple, Union

from swanlab.exceptions import DataStoreError
from swanlab.sdk.internal.core_python.store import (
    LEVELDBLOG_BLOCK_LEN,
    LEVELDBLOG_FIRST,
    LEVELDBLOG_FULL,
    LEVELDBLOG_HEADER_LEN,
    LEVELDBLOG_LAST,
    LEVELDBLOG_MIDDLE,
)

__all__ = ["WBDataStoreReader", "DataStoreError", "DataStoreCorruptionError", "DataStoreTruncatedError"]

LEVELDBLOG_TYPES = (LEVELDBLOG_FULL, LEVELDBLOG_FIRST, LEVELDBLOG_MIDDLE, LEVELDBLOG_LAST)

# W&B 事务日志文件头：ident(4) + magic(uint16 LE) + version(uint8)
WANDB_HEADER_IDENT = b":W&B"
WANDB_HEADER_MAGIC = 0xBEE1  # zlib.crc32(b"Weights & Biases") & 0xffff
WANDB_HEADER_VERSION = 0
# 预计算各分片类型字节的 crc32，作为数据校验和的种子（键即类型值）
_CRC = {t: zlib.crc32(bytes([t])) & 0xFFFFFFFF for t in LEVELDBLOG_TYPES}


class DataStoreCorruptionError(DataStoreError):
    """文件内容违反事务日志格式（坏文件头、校验和不匹配、非法分片类型或填充）。"""


class DataStoreTruncatedError(DataStoreError):
    """文件在记录中途结束（可能仍在写入，或上次写入被中断）。"""


class WBDataStoreReader:
    """W&B leveldb 事务日志只读扫描器。独立实现，不继承 DataStoreReader。"""

    def __init__(self) -> None:
        self._fp: Optional[IO[bytes]] = None
        self._index: int = 0

    def __enter__(self) -> "WBDataStoreReader":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def open(self, filename: Union[Path, str]) -> None:
        """打开文件并校验 W&B 文件头；校验失败时关闭文件句柄。"""
        if self._fp is not None:
            raise DataStoreError("reader is already open")
        self._fp = open(filename, "rb")
        self._index = 0
        try:
            self._read_header()
        except Exception:
            self._fp.close()
            self._fp = None
            raise

    def close(self) -> None:
        """幂等关闭；open 失败后调用也安全。"""
        if self._fp is not None:
            self._fp.close()
            self._fp = None

    def scan(self) -> Optional[bytes]:
        """返回下一条完整记录；文件结束时返回 None。跨块分片自动重组。"""
        dtype, data = self._read_record_strict()
        if dtype is None:
            return None
        if dtype == LEVELDBLOG_FULL:
            return data
        if dtype != LEVELDBLOG_FIRST:
            raise DataStoreCorruptionError(f"chunk type {dtype} ending at offset {self._index}: expected FULL or FIRST")

        chunks = [data]
        while True:
            dtype, chunk = self._read_record_strict()
            if dtype is None:
                raise DataStoreTruncatedError(f"fragmented record ends at offset {self._index} without a LAST chunk")
            chunks.append(chunk)
            if dtype == LEVELDBLOG_LAST:
                return b"".join(chunks)
            if dtype != LEVELDBLOG_MIDDLE:
                raise DataStoreCorruptionError(
                    f"chunk type {dtype} ending at offset {self._index}: expected MIDDLE or LAST"
                )

    def _read_header(self) -> None:
        assert self._fp is not None
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        if len(header) != LEVELDBLOG_HEADER_LEN:
            raise DataStoreTruncatedError(f"truncated file header: expected {LEVELDBLOG_HEADER_LEN} bytes")
        ident, magic, version = struct.unpack("<4sHB", header)
        if ident != WANDB_HEADER_IDENT or magic != WANDB_HEADER_MAGIC:
            raise DataStoreCorruptionError(f"invalid W&B file header: ident={ident!r}, magic={magic:#06x}")
        if version != WANDB_HEADER_VERSION:
            raise DataStoreError(
                f"unsupported W&B transaction log version {version} "
                f"(supported: {WANDB_HEADER_VERSION}); please upgrade SwanLab to read it"
            )
        self._index += len(header)

    def _read_record_strict(self) -> Tuple[Optional[int], bytes]:
        """读取单个分片并校验分帧；文件干净结束时返回 (None, b"")。"""
        if self._fp is None:
            raise DataStoreError("file not open for scanning; call open() first")
        self._skip_block_trailer()

        offset = self._index
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        if len(header) == 0:
            return None, b""
        if len(header) != LEVELDBLOG_HEADER_LEN:
            raise DataStoreTruncatedError(f"truncated chunk header at offset {offset}")

        checksum, length, dtype = struct.unpack("<IHB", header)
        if dtype not in LEVELDBLOG_TYPES:
            raise DataStoreCorruptionError(f"invalid chunk type {dtype} at offset {offset}")
        # 分片不得跨越块边界：超界长度按损坏处理，避免级联读入垃圾
        room_in_block = LEVELDBLOG_BLOCK_LEN - LEVELDBLOG_HEADER_LEN - (offset % LEVELDBLOG_BLOCK_LEN)
        if length > room_in_block:
            raise DataStoreCorruptionError(f"chunk length {length} crosses block boundary at offset {offset}")

        data = self._fp.read(length)
        if len(data) != length:
            raise DataStoreTruncatedError(f"truncated chunk payload at offset {offset}: {len(data)}/{length} bytes")
        if zlib.crc32(data, _CRC[dtype]) & 0xFFFFFFFF != checksum:
            raise DataStoreCorruptionError(f"chunk checksum mismatch at offset {offset}")

        self._index = offset + LEVELDBLOG_HEADER_LEN + length
        return dtype, data

    def _skip_block_trailer(self) -> None:
        """跳过块尾不足一个分片头的零填充。"""
        assert self._fp is not None  # 仅供 _read_record_strict 在守卫后调用
        space_left = LEVELDBLOG_BLOCK_LEN - (self._index % LEVELDBLOG_BLOCK_LEN)
        if space_left >= LEVELDBLOG_HEADER_LEN:
            return
        offset = self._index
        pad = self._fp.read(space_left)
        if pad == b"\x00" * space_left:
            self._index += space_left  # 完整填充，进入下一块
        elif not pad:
            return  # 文件止于记录边界
        elif pad == b"\x00" * len(pad):
            raise DataStoreTruncatedError(f"incomplete block padding at offset {offset}")
        else:
            raise DataStoreCorruptionError(f"invalid block padding at offset {offset}: {pad!r}")
