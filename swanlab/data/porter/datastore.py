"""
@author: cunyue
@file: datastore.py
@time: 2025/6/5 16:32
@description: 记录的数据遵循LevelDB格式：https://github.com/google/leveldb/blob/main/doc/log_format.md
我们使用 crc32 计算数据校验和，crc32 相对轻量，且计算速度较快
字符编码使用 utf-8, 确保数据兼容性 （这与 LevelDB 的规范有所冲突，如果有必要，未来可以升级版本并通过LEVELDBLOG_HEADER_VERSION兼容）
这为后续引入 protobuf 或其他序列化格式打下基础
DataStore 大致代码借鉴自 W&B

文件头版本区分了不同更新时的格式变化，这方面我们不会做向下兼容，即低版本文件不一定能被新版本读取，可以通过版本降级来区分不同版本号
"""

import os
import struct
import zlib
from typing import Optional, Any, IO, Tuple

from swanlab.error import ValidationError

LEVELDBLOG_HEADER_LEN = 7
LEVELDBLOG_BLOCK_LEN = 32768
LEVELDBLOG_DATA_LEN = LEVELDBLOG_BLOCK_LEN - LEVELDBLOG_HEADER_LEN

LEVELDBLOG_FULL = 1
LEVELDBLOG_FIRST = 2
LEVELDBLOG_MIDDLE = 3
LEVELDBLOG_LAST = 4


LEVELDBLOG_HEADER_IDENT = ":SWL"
LEVELDBLOG_HEADER_MAGIC = 0xE1D6  # zlib.crc32(bytes("SwanLab", 'utf-8')) & 0xffff
LEVELDBLOG_HEADER_VERSION = 1


def strtobytes(x):
    """
    文件转字符串
    """
    return bytes(x, "utf-8")


def bytestostr(x):
    return str(x, 'utf-8')


class DataStore:

    def __init__(self):
        self._filename: Optional[str] = None
        self._fp: Optional[IO[Any]] = None
        # 当前文件的偏移量
        self._index: int = 0
        # 当前文件的已刷写偏移量
        self._flush_offset = 0
        # 日志系统预计算并缓存CRC32校验值，缓存每一个数据类型的CRC32值，分别存在各自的索引位置
        self._crc = [0] * (LEVELDBLOG_LAST + 1)
        for x in range(1, LEVELDBLOG_LAST + 1):
            self._crc[x] = zlib.crc32(strtobytes(chr(x))) & 0xFFFFFFFF

        # 是否为扫描模式打开文件
        self._opened_for_scan = False
        # 当前文件大小（仅在扫描模式下有效）
        self._size_bytes: int = 0

    # ---------------------------------- 读取 ----------------------------------

    def open_for_scan(self, filename: str):
        self._filename = filename
        self._fp = open(filename, "r+b")
        self._index = 0
        self._size_bytes = os.stat(filename).st_size
        self._opened_for_scan = True
        self._read_header()

    def _read_header(self):
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        assert (
            len(header) == LEVELDBLOG_HEADER_LEN
        ), f"header is {len(header)} bytes instead of the expected {LEVELDBLOG_HEADER_LEN}"
        ident, magic, version = struct.unpack("<4sHB", header)
        if ident != strtobytes(LEVELDBLOG_HEADER_IDENT):
            raise Exception("Invalid header")
        if magic != LEVELDBLOG_HEADER_MAGIC:
            raise Exception("Invalid header")
        if version != LEVELDBLOG_HEADER_VERSION:
            raise Exception(
                f"Invalid backup version: {version}. For supported versions, see: https://docs.swanlab.cn/api/cli-swanlab-sync.html"
            )
        self._index += len(header)

    def _scan_record(self) -> Optional[Tuple[int, bytes]]:
        """
        扫描一条记录
        """
        assert self._opened_for_scan, "file not open for scanning"
        # 1. 读取数据头
        header = self._fp.read(LEVELDBLOG_HEADER_LEN)
        if len(header) == 0:
            return None
        assert (
            len(header) == LEVELDBLOG_HEADER_LEN
        ), f"record header is {len(header)} bytes instead of the expected {LEVELDBLOG_HEADER_LEN}"
        # 2. 解析数据头并校验数据完整性
        checksum, data_length, data_type = struct.unpack("<IHB", header)
        self._index += LEVELDBLOG_HEADER_LEN
        data = self._fp.read(data_length)
        checksum_computed = zlib.crc32(data, self._crc[data_type]) & 0xFFFFFFFF
        if checksum != checksum_computed:
            raise ValidationError("Invalid record checksum, data may be corrupt")
        self._index += data_length
        # 3. 返回数据
        return int(data_type), data

    def scan(self) -> Optional[str]:
        """
        扫描日志文件，返回一条记录
        """
        # 1. 一次读取一条记录，如果剩余空间不足存储数据头，校验并跳过，此为写入的逆操作
        offset = self._index % LEVELDBLOG_BLOCK_LEN
        space_left = LEVELDBLOG_BLOCK_LEN - offset
        if space_left < LEVELDBLOG_HEADER_LEN:
            pad_check = strtobytes("\x00" * space_left)
            pad = self._fp.read(space_left)
            # 校验必须为0
            assert pad == pad_check, "invalid padding"
            self._index += space_left
        # 2. 扫描一条记录
        record = self._scan_record()
        if record is None:  # eof
            return None
        dtype, data = record
        if dtype == LEVELDBLOG_FULL:
            return bytestostr(data)
        # 3. 如果是第一条记录，则继续扫描直到找到最后一条记录
        assert dtype == LEVELDBLOG_FIRST, f"expected record to be type {LEVELDBLOG_FIRST} but found {dtype}"
        while True:
            record = self._scan_record()
            if record is None:  # eof
                return None
            dtype, new_data = record
            if dtype == LEVELDBLOG_LAST:
                data += new_data
                break
            assert dtype == LEVELDBLOG_MIDDLE, f"expected record to be type {LEVELDBLOG_MIDDLE} but found {dtype}"
            data += new_data
        return bytestostr(data)

    def __iter__(self):
        """
        实现迭代器接口，允许使用 for 循环遍历日志文件，仅在文件已打开并且处于扫描模式时有效
        """
        assert self._opened_for_scan, "file not open for scanning, cannot iterate"
        return self

    def __next__(self):
        record = self.scan()
        if record is None:
            raise StopIteration("End of file reached")
        return record

    # ---------------------------------- 写入 ----------------------------------

    def open_for_write(self, filename: str):
        self._filename = filename
        self._fp = open(filename, "xb")
        # 写入文件头, 长度等于 LEVELDBLOG_HEADER_LEN
        data = struct.pack(
            "<4sHB",
            strtobytes(LEVELDBLOG_HEADER_IDENT),
            LEVELDBLOG_HEADER_MAGIC,
            LEVELDBLOG_HEADER_VERSION,
        )
        assert len(data) == LEVELDBLOG_HEADER_LEN, f"header size is {len(data)} bytes, expected {LEVELDBLOG_HEADER_LEN}"
        self._fp.write(data)
        self._index += len(data)

    def _write_record(self, data: bytes, data_type: int = LEVELDBLOG_FULL):
        """
        写入记录到日志文件
        """
        assert len(data) + LEVELDBLOG_HEADER_LEN <= (
            LEVELDBLOG_BLOCK_LEN - self._index % LEVELDBLOG_BLOCK_LEN
        ), "not enough space to write new records"
        data_length = len(data)
        # 计算校验值，校验值为对数据和数据类型的 CRC32 校验和
        checksum = zlib.crc32(data, self._crc[data_type]) & 0xFFFFFFFF
        # 写入数据头，格式为：<IHB>，分别表示校验和、数据长度和数据类型
        # I: unsigned int (4 bytes), H: unsigned short (2 bytes), B: unsigned char (1 byte)
        self._fp.write(struct.pack("<IHB", checksum, data_length, data_type))
        if data_length:
            self._fp.write(data)
        self._index += LEVELDBLOG_HEADER_LEN + len(data)

    def write(self, s: str):
        """
        写入数据到日志文件，遵循 LevelDB 规范
        :param s: 要写入的数据，必须是字符串形式
        :return: 返回写入的起始偏移量、当前偏移量和已刷写偏移量
        """
        data = strtobytes(s)
        # 1. 计算偏移量
        start_offset = self._index
        offset = self._index % LEVELDBLOG_BLOCK_LEN
        space_left = LEVELDBLOG_BLOCK_LEN - offset
        data_used = 0
        data_left = len(data)
        # 2. 剩余长度小于数据头长度则填充0，归位到下一个块
        if space_left < LEVELDBLOG_HEADER_LEN:
            pad = "\x00" * space_left
            self._fp.write(strtobytes(pad))
            self._index += space_left
            space_left = LEVELDBLOG_BLOCK_LEN
        # 3. 如果剩余长度大于等于数据长度，则直接写入
        if data_left + LEVELDBLOG_HEADER_LEN <= space_left:
            self._write_record(data)
        # 4. 否则需要分块写入（注意此时我们可能在一个块的中间）
        else:
            # 4.1 写入第一个数据块，确保接下来数据独占一个块
            data_room = space_left - LEVELDBLOG_HEADER_LEN
            self._write_record(data[:data_room], LEVELDBLOG_FIRST)
            data_used += data_room
            data_left -= data_room
            assert data_left, "data_left should be non-zero"
            # 4.2 写入中间数据
            while data_left > LEVELDBLOG_DATA_LEN:
                self._write_record(
                    data[data_used : data_used + LEVELDBLOG_DATA_LEN],
                    LEVELDBLOG_MIDDLE,
                )
                data_used += LEVELDBLOG_DATA_LEN
                data_left -= LEVELDBLOG_DATA_LEN
            # 4.3 写入最后一个数据块
            self._write_record(data[data_used:], LEVELDBLOG_LAST)
            # 刷写完整数据
            try:
                self._fp.flush()
                os.fsync(self._fp.fileno())
            except OSError:
                # 如果操作系统不支持 fsync，可能会抛出 OSError，忽略此错误即可
                pass
            self._flush_offset = self._index

        return start_offset, self._index, self._flush_offset

    # ---------------------------------- 辅助函数 ----------------------------------

    def ensure_flushed(self) -> None:
        self._fp.flush()

    def close(self):
        # 关闭文件句柄
        self._fp.close()
