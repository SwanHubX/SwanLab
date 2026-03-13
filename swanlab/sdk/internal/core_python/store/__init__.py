"""
@author: cunyue
@file: __init__.py
@time: 2025/6/5 16:32
@description: 记录的数据遵循LevelDB格式：https://github.com/google/leveldb/blob/main/doc/log_format.md
我们使用 crc32 计算数据校验和，crc32 相对轻量，且计算速度较快
DataStore 大致代码借鉴自 W&B

文件头版本区分了不同更新时的格式变化，这方面我们不会做向下兼容，即低版本文件不一定能被新版本读取，可以通过版本降级来区分不同版本号
"""

import os
import struct
import zlib
from typing import IO, Any, Optional, Tuple

from swanlab.exceptions import DataStoreError

__all__ = ["DataStoreWriter", "DataStoreReader", "DataStoreError"]

LEVELDBLOG_HEADER_LEN = 7
LEVELDBLOG_BLOCK_LEN = 32768
LEVELDBLOG_DATA_LEN = LEVELDBLOG_BLOCK_LEN - LEVELDBLOG_HEADER_LEN

LEVELDBLOG_FULL = 1
LEVELDBLOG_FIRST = 2
LEVELDBLOG_MIDDLE = 3
LEVELDBLOG_LAST = 4

LEVELDBLOG_HEADER_IDENT = b":SWL"
LEVELDBLOG_HEADER_MAGIC = 0xE1D6  # zlib.crc32(bytes("SwanLab", 'utf-8')) & 0xffff
LEVELDBLOG_HEADER_VERSION = 1

# 模块级 CRC 预计算，构造一次，读写两侧共用
_CRC = [0] * (LEVELDBLOG_LAST + 1)
for _x in range(1, LEVELDBLOG_LAST + 1):
    _CRC[_x] = zlib.crc32(bytes([_x])) & 0xFFFFFFFF


# ===========================================================================
# 写入
# ===========================================================================


class DataStoreWriter:
    """追加写入器，持有一个长期打开的二进制文件句柄。"""

    def __init__(self):
        self._fp: Optional[IO[Any]] = None
        self._index: int = 0
        self._flush_offset: int = 0

    def open(self, filename: str) -> None:
        """创建并初始化文件，文件已存在时抛出 FileExistsError。"""
        self._fp = open(filename, "xb")
        header = struct.pack("<4sHB", LEVELDBLOG_HEADER_IDENT, LEVELDBLOG_HEADER_MAGIC, LEVELDBLOG_HEADER_VERSION)
        assert len(header) == LEVELDBLOG_HEADER_LEN
        self._fp.write(header)
        self._index += len(header)

    def write(self, data: bytes) -> None:
        """写入任意字节，遵循 LevelDB log 分块规范。"""
        assert self._fp is not None, "writer is not open"
        offset = self._index % LEVELDBLOG_BLOCK_LEN
        space_left = LEVELDBLOG_BLOCK_LEN - offset
        data_used = 0
        data_left = len(data)
        # 剩余空间不足一个 header：填充 0，归位到下一个块
        if space_left < LEVELDBLOG_HEADER_LEN:
            self._fp.write(b"\x00" * space_left)
            self._index += space_left
            space_left = LEVELDBLOG_BLOCK_LEN
        # 数据可以放入当前块
        if data_left + LEVELDBLOG_HEADER_LEN <= space_left:
            self._write_record(data)
        # 否则分块写入
        else:
            data_room = space_left - LEVELDBLOG_HEADER_LEN
            self._write_record(data[:data_room], LEVELDBLOG_FIRST)
            data_used += data_room
            data_left -= data_room
            assert data_left, "data_left should be non-zero"
            while data_left > LEVELDBLOG_DATA_LEN:
                self._write_record(data[data_used : data_used + LEVELDBLOG_DATA_LEN], LEVELDBLOG_MIDDLE)
                data_used += LEVELDBLOG_DATA_LEN
                data_left -= LEVELDBLOG_DATA_LEN
            self._write_record(data[data_used:], LEVELDBLOG_LAST)
        # 每次 write 后统一 fsync，保证落盘
        try:
            self._fp.flush()
            os.fsync(self._fp.fileno())
        except OSError:
            pass
        self._flush_offset = self._index

    def ensure_flushed(self) -> None:
        assert self._fp is not None, "writer is not open"
        self._fp.flush()

    def close(self) -> None:
        assert self._fp is not None, "writer is not open"
        self._fp.flush()
        self._fp.close()
        self._fp = None

    def _write_record(self, data: bytes, data_type: int = LEVELDBLOG_FULL) -> None:
        assert self._fp is not None
        assert len(data) + LEVELDBLOG_HEADER_LEN <= (LEVELDBLOG_BLOCK_LEN - self._index % LEVELDBLOG_BLOCK_LEN), (
            "not enough space to write new records"
        )
        checksum = zlib.crc32(data, _CRC[data_type]) & 0xFFFFFFFF
        self._fp.write(struct.pack("<IHB", checksum, len(data), data_type))
        if data:
            self._fp.write(data)
        self._index += LEVELDBLOG_HEADER_LEN + len(data)


# ===========================================================================
# 读取
# ===========================================================================


class DataStoreReader:
    """顺序扫描读取器，实现迭代器协议。"""

    def __init__(self):
        self._fp: Optional[IO[Any]] = None
        self._index: int = 0

    def open(self, filename: str) -> None:
        """打开文件并校验文件头。"""
        self._fp = open(filename, "rb")
        self._index = 0
        self._read_header()

    def scan(self) -> Optional[bytes]:
        """读取下一条完整记录，到达文件末尾时返回 None。"""
        assert self._fp is not None, "reader is not open"
        # 剩余空间不足一个 header：跳过填充字节
        offset = self._index % LEVELDBLOG_BLOCK_LEN
        space_left = LEVELDBLOG_BLOCK_LEN - offset
        if space_left < LEVELDBLOG_HEADER_LEN:
            pad = self._fp.read(space_left)
            assert pad == b"\x00" * space_left, "invalid padding"
            self._index += space_left
        record = self._read_record()
        if record is None:
            return None
        dtype, data = record
        if dtype == LEVELDBLOG_FULL:
            return data
        # 分块记录：连续读取直到 LAST
        assert dtype == LEVELDBLOG_FIRST, f"expected record type {LEVELDBLOG_FIRST}, got {dtype}"
        while True:
            record = self._read_record()
            if record is None:
                return None
            dtype, chunk = record
            data += chunk
            if dtype == LEVELDBLOG_LAST:
                return data
            assert dtype == LEVELDBLOG_MIDDLE, f"expected record type {LEVELDBLOG_MIDDLE}, got {dtype}"

    def close(self) -> None:
        assert self._fp is not None, "reader is not open"
        self._fp.close()
        self._fp = None

    def __iter__(self):
        assert self._fp is not None, "reader is not open"
        return self

    def __next__(self) -> bytes:
        record = self.scan()
        if record is None:
            raise StopIteration
        return record

    def _read_header(self) -> None:
        assert self._fp is not None
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        assert len(header) == LEVELDBLOG_HEADER_LEN, (
            f"header is {len(header)} bytes instead of the expected {LEVELDBLOG_HEADER_LEN}"
        )
        ident, magic, version = struct.unpack("<4sHB", header)
        if ident != LEVELDBLOG_HEADER_IDENT:
            raise DataStoreError("Invalid header ident")
        if magic != LEVELDBLOG_HEADER_MAGIC:
            raise DataStoreError("Invalid header magic")
        if version != LEVELDBLOG_HEADER_VERSION:
            raise DataStoreError(
                f"Invalid run version: {version}. For supported versions, see: https://docs.swanlab.cn/api/cli-swanlab-sync.html"
            )
        self._index += len(header)

    def _read_record(self) -> Optional[Tuple[int, bytes]]:
        assert self._fp is not None
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        if len(header) == 0:
            return None
        assert len(header) == LEVELDBLOG_HEADER_LEN, (
            f"record header is {len(header)} bytes instead of the expected {LEVELDBLOG_HEADER_LEN}"
        )
        checksum, data_length, data_type = struct.unpack("<IHB", header)
        self._index += LEVELDBLOG_HEADER_LEN
        data = self._fp.read(data_length)
        if zlib.crc32(data, _CRC[data_type]) & 0xFFFFFFFF != checksum:
            raise DataStoreError("Invalid record checksum, data may be corrupt")
        self._index += data_length
        return int(data_type), data
