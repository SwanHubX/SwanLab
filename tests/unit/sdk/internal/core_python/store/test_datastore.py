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

    def test_legacy_version_rejected(self, tmp_path: Path):
        """旧 SDK（v1）生成的文件应被拒绝，并提示版本不匹配。"""
        p = tmp_path / "legacy.swanlab"
        _write_raw_header(p, LEVELDBLOG_HEADER_IDENT, LEVELDBLOG_HEADER_MAGIC, 1)
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

    def test_iteration_stops_at_corrupted_record(self, tmp_path: Path):
        p = tmp_path / "test.swanlab"
        first = b"first"
        second = b"second"
        third = b"third"
        write_records(p, first, second, third)

        second_data_offset = LEVELDBLOG_HEADER_LEN + LEVELDBLOG_HEADER_LEN + len(first) + LEVELDBLOG_HEADER_LEN
        data = bytearray(p.read_bytes())
        data[second_data_offset] ^= 0xFF
        p.write_bytes(bytes(data))

        r = DataStoreReader()
        r.open(str(p))
        assert list(r) == [first]
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


# ---------------------------------------------------------------------------
# write / sync 分离：fsync 频率与落盘语义
# ---------------------------------------------------------------------------


class TestWriteSyncSeparation:
    """write() 仅写缓冲区不 fsync，sync() 才落盘；批量写后一次 sync。"""

    def test_write_does_not_fsync(self, tmp_path: Path, mocker):
        """单纯 write 不应触发 os.fsync。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        w.write(b"a")
        w.write(b"b")
        assert spy.call_count == 0
        w.close()

    def test_sync_triggers_single_fsync(self, tmp_path: Path, mocker):
        """连续多次 write 后，一次 sync 只产生一次 fsync。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        for i in range(50):
            w.write(f"record_{i}".encode())
        assert spy.call_count == 0
        w.sync()
        assert spy.call_count == 1
        w.close()

    def test_sync_without_dirty_is_noop(self, tmp_path: Path, mocker):
        """无脏数据时 sync 不应调用 fsync（含重复 sync）。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        # 刚 open、未 write：sync 是 no-op
        w.sync()
        assert spy.call_count == 0
        # write 后 sync 一次，再次 sync 不重复 fsync
        w.write(b"x")
        w.sync()
        w.sync()
        assert spy.call_count == 1
        w.close()

    def test_close_flushes_pending_writes(self, tmp_path: Path, mocker):
        """write 后未显式 sync，close 也应保证落盘（close 内含 sync）。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        w.write(b"only via close")
        assert spy.call_count == 0
        w.close()
        assert spy.call_count == 1
        # 真实落盘内容可读回
        assert read_all(p) == [b"only via close"]

    def test_close_without_dirty_does_not_fsync(self, tmp_path: Path, mocker):
        """已 sync 过、无新写入时 close 不应再 fsync。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        w.write(b"data")
        w.sync()
        assert spy.call_count == 1
        w.close()
        assert spy.call_count == 1

    def test_ensure_flushed_delegates_to_sync(self, tmp_path: Path, mocker):
        """ensure_flushed 等价于 sync，会落盘脏数据。"""
        spy = mocker.patch("swanlab.sdk.internal.core_python.store.os.fsync")
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        w.write(b"data")
        w.ensure_flushed()
        assert spy.call_count == 1
        w.close()

    def test_data_in_buffer_not_visible_before_sync(self, tmp_path: Path):
        """write/open 后未 sync 时，数据（连文件头）都停留在 Python 用户态缓冲区。

        新开一个独立 reader 句柄此时连文件头都读不到；sync 后才推给 OS 可读回。
        写入仅入用户态 buffer，sync() 才落盘。
        """
        p = tmp_path / "test.swanlab"
        w = DataStoreWriter()
        w.open(str(p))
        w.write(b"buffered record")
        # 此时另开句柄读取：连文件头都还在缓冲区，reader 校验文件头失败
        r = DataStoreReader()
        with pytest.raises(AssertionError):
            r.open(str(p))
        # sync 后数据推给 OS，可被完整读回
        w.sync()
        assert read_all(p) == [b"buffered record"]
        w.close()

    def test_batch_write_roundtrip(self, tmp_path: Path):
        """批量 write + 一次 sync 后，全部记录可按序读回。"""
        p = tmp_path / "test.swanlab"
        payloads = [f"r{i}".encode() for i in range(200)]
        w = DataStoreWriter()
        w.open(str(p))
        for payload in payloads:
            w.write(payload)
        w.sync()
        w.close()
        assert read_all(p) == payloads
