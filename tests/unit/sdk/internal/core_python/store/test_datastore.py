"""
@author: cunyue
@file: test_datastore.py
@time: 2026/3/13
@description: 测试 DataStoreWriter / DataStoreReader（LevelDB log 格式）的读写行为
"""

import struct
from pathlib import Path

import pytest

from swanlab.sdk.internal.core_python.store import (
    LEVELDBLOG_BLOCK_LEN,
    LEVELDBLOG_DATA_LEN,
    LEVELDBLOG_HEADER_IDENT,
    LEVELDBLOG_HEADER_LEN,
    LEVELDBLOG_HEADER_MAGIC,
    LEVELDBLOG_HEADER_VERSION,
    DataStoreError,
    DataStoreReader,
    DataStoreWriter,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def write_records(path: Path, *payloads: bytes) -> None:
    w = DataStoreWriter()
    w.open(str(path))
    for payload in payloads:
        w.write(payload)
    w.close()


def read_all(path: Path) -> list[bytes]:
    r = DataStoreReader()
    r.open(str(path))
    results = list(r)
    r.close()
    return results


# ---------------------------------------------------------------------------
# 写入 / 读取基本行为
# ---------------------------------------------------------------------------


class TestWriteAndScan:
    def test_single_record(self, tmp_path: Path):
        p = tmp_path / "test.swanlab"
        payload = b"hello swanlab"
        write_records(p, payload)
        assert read_all(p) == [payload]

    def test_multiple_records(self, tmp_path: Path):
        p = tmp_path / "test.swanlab"
        payloads = [b"record_1", b"record_2", b"record_3"]
        write_records(p, *payloads)
        assert read_all(p) == payloads

    def test_binary_payload(self, tmp_path: Path):
        """任意二进制（如 protobuf 序列化结果）可以原样读回。"""
        p = tmp_path / "test.swanlab"
        payload = bytes(range(256))
        write_records(p, payload)
        assert read_all(p) == [payload]

    def test_scan_returns_none_at_eof(self, tmp_path: Path):
        p = tmp_path / "test.swanlab"
        write_records(p, b"only_one")
        r = DataStoreReader()
        r.open(str(p))
        r.scan()  # 消费唯一一条
        assert r.scan() is None
        r.close()

    def test_empty_file_scan_returns_none(self, tmp_path: Path):
        """只有文件头、无数据记录时，scan 应立即返回 None。"""
        p = tmp_path / "test.swanlab"
        write_records(p)
        assert read_all(p) == []


# ---------------------------------------------------------------------------
# 大数据 / 跨块写入
# ---------------------------------------------------------------------------


class TestMultiBlockWrite:
    def test_record_spanning_two_blocks(self, tmp_path: Path):
        """数据大于单块可用空间，触发 FIRST / LAST 分块写入。"""
        p = tmp_path / "test.swanlab"
        payload = b"x" * (LEVELDBLOG_DATA_LEN + 1)
        write_records(p, payload)
        assert read_all(p) == [payload]

    def test_record_spanning_three_blocks(self, tmp_path: Path):
        """超过两块，触发 FIRST / MIDDLE / LAST 路径。"""
        p = tmp_path / "test.swanlab"
        payload = b"y" * (LEVELDBLOG_DATA_LEN * 2 + 1)
        write_records(p, payload)
        assert read_all(p) == [payload]

    def test_block_boundary_padding(self, tmp_path: Path):
        """写入数据使剩余空间不足一个 header，应插入 0 填充后继续写下一条。"""
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        # 先写一条撑到距块末尾恰好 < HEADER_LEN 的位置
        # 文件头占 7 字节，故首条 payload 大小 = BLOCK_LEN - HEADER_LEN*2 - 7
        first_len = LEVELDBLOG_BLOCK_LEN - LEVELDBLOG_HEADER_LEN * 2 - LEVELDBLOG_HEADER_LEN
        w.write(b"a" * first_len)
        second = b"b" * 10
        w.write(second)
        w.close()

        results = read_all(p)
        assert results[0] == b"a" * first_len
        assert results[1] == second


# ---------------------------------------------------------------------------
# 迭代器协议
# ---------------------------------------------------------------------------


class TestIterator:
    def test_for_loop(self, tmp_path: Path):
        p = tmp_path / "test.swanlab"
        payloads = [b"alpha", b"beta", b"gamma"]
        write_records(p, *payloads)
        r = DataStoreReader()
        r.open(str(p))
        result = list(r)
        r.close()
        assert result == payloads

    def test_iter_before_open_raises(self):
        """未 open 时，iter 应抛出 AssertionError。"""
        r = DataStoreReader()
        with pytest.raises(AssertionError):
            iter(r)


# ---------------------------------------------------------------------------
# 文件头校验
# ---------------------------------------------------------------------------


def _write_raw_header(path: Path, ident: bytes, magic: int, version: int):
    with open(path, "wb") as f:
        f.write(struct.pack("<4sHB", ident, magic, version))


class TestHeaderValidation:
    def test_invalid_ident(self, tmp_path: Path):
        p = tmp_path / "bad.swanlab"
        _write_raw_header(p, b"XXXX", LEVELDBLOG_HEADER_MAGIC, LEVELDBLOG_HEADER_VERSION)
        r = DataStoreReader()
        with pytest.raises(DataStoreError, match="ident"):
            r.open(str(p))

    def test_invalid_magic(self, tmp_path: Path):
        p = tmp_path / "bad.swanlab"
        _write_raw_header(p, LEVELDBLOG_HEADER_IDENT, 0x0000, LEVELDBLOG_HEADER_VERSION)  # type: ignore
        r = DataStoreReader()
        with pytest.raises(DataStoreError, match="magic"):
            r.open(str(p))

    def test_invalid_version(self, tmp_path: Path):
        p = tmp_path / "bad.swanlab"
        _write_raw_header(p, LEVELDBLOG_HEADER_IDENT, LEVELDBLOG_HEADER_MAGIC, 99)  # type: ignore
        r = DataStoreReader()
        with pytest.raises(DataStoreError, match="version"):
            r.open(str(p))


# ---------------------------------------------------------------------------
# CRC 校验
# ---------------------------------------------------------------------------


class TestChecksumValidation:
    def test_corrupted_data_raises(self, tmp_path: Path):
        """翻转数据区某个字节，scan 时应抛出 DataStoreError。"""
        p = tmp_path / "test.swanlab"
        write_records(p, b"intact data")
        # 文件头 + record header 之后的第一个数据字节处注入错误
        corrupt_offset = LEVELDBLOG_HEADER_LEN * 2
        data = bytearray(p.read_bytes())
        data[corrupt_offset] ^= 0xFF
        p.write_bytes(bytes(data))

        r = DataStoreReader()
        r.open(str(p))
        with pytest.raises(DataStoreError, match="checksum"):
            r.scan()
        r.close()


# ---------------------------------------------------------------------------
# 文件管理
# ---------------------------------------------------------------------------


class TestFileHandling:
    def test_open_for_write_exclusive(self, tmp_path: Path):
        """文件已存在时，open 应抛出 FileExistsError。"""
        p = tmp_path / "test.swanlab"
        p.touch()
        w = DataStoreWriter()
        with pytest.raises(FileExistsError):
            w.open(str(p))

    def test_data_readable_after_close(self, tmp_path: Path):
        """close() 后数据应已落盘，可被 reader 完整读取。"""
        p = tmp_path / "test.swanlab"
        payload = b"persist me"
        write_records(p, payload)
        assert read_all(p) == [payload]

    def test_write_before_open_raises(self):
        """未 open 直接 write 应抛出 AssertionError。"""
        w = DataStoreWriter()
        with pytest.raises(AssertionError):
            w.write(b"data")

    def test_scan_before_open_raises(self):
        """未 open 直接 scan 应抛出 AssertionError。"""
        r = DataStoreReader()
        with pytest.raises(AssertionError):
            r.scan()
