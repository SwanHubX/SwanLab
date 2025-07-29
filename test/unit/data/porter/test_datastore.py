"""
@author: cunyue
@file: test_datastore.py
@time: 2025/6/5 18:33
@description: 测试 DataStore 类的功能，包括数据写入、读取、校验和计算等
"""

import struct

import pytest
from nanoid import generate

from swanlab.data.porter.datastore import DataStore, LEVELDBLOG_HEADER_LEN, LEVELDBLOG_BLOCK_LEN, LEVELDBLOG_MIDDLE
from swanlab.error import ValidationError

logs = [generate(size=l) for l in range(1, 100001, 1000)]


@pytest.fixture
def temp_log_file(tmp_path):
    """创建临时日志文件并写入测试数据"""
    log_file = tmp_path / "test.log"
    ds = DataStore()
    ds.open_for_write(str(log_file))
    ds.write("Valid record 1")
    ds.write("Another valid record")
    ds.close()
    return log_file


def test_write(temp_log_file):
    """
    测试文件写入功能
    """
    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    records = [next(ds), next(ds)]
    assert records == ["Valid record 1", "Another valid record"]


def test_scan(temp_log_file):
    """
    测试文件扫描功能
    """
    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    records = [next(ds), next(ds)]
    assert records == ["Valid record 1", "Another valid record"]


def test_tampered_data_content(temp_log_file):
    """篡改数据内容应触发校验失败"""
    # 篡改第一条记录的数据部分
    with open(temp_log_file, "r+b") as f:
        # 定位到第一条记录的数据起始位置（跳过文件头）
        f.seek(LEVELDBLOG_HEADER_LEN + 7)  # 7=CRC(4)+Length(2)+Type(1)
        _ = f.read(13)  # "Valid record 1"的长度
        f.seek(LEVELDBLOG_HEADER_LEN + 7)
        f.write(b"Tampered!!!!")  # 篡改数据

    # 尝试读取应触发异常
    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    with pytest.raises(ValidationError, match="Invalid record checksum"):
        next(ds)  # 读取第一条记录


def test_tampered_checksum(temp_log_file):
    """篡改校验和字段应触发校验失败"""
    with open(temp_log_file, "r+b") as f:
        # 定位到第一条记录的CRC字段（跳过文件头）
        f.seek(LEVELDBLOG_HEADER_LEN)
        f.write(struct.pack("<I", 0xDEADBEEF))  # 写入无效CRC

    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    with pytest.raises(ValidationError, match="Invalid record checksum"):
        next(ds)


def test_tampered_record_type(temp_log_file):
    """篡改记录类型应触发校验失败"""
    with open(temp_log_file, "r+b") as f:
        # 定位到记录类型字段（CRC+Length之后）
        f.seek(LEVELDBLOG_HEADER_LEN + 6)  # 4(CRC)+2(Length)
        f.write(bytes([LEVELDBLOG_MIDDLE]))  # 改为无效类型

    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    with pytest.raises(ValidationError, match="Invalid record checksum"):
        next(ds)


def test_tampered_multi_block_record(tmp_path):
    """测试多块记录的篡改场景"""
    # 创建临时日志文件
    # 不使用 pytest 自带的 tmp_log_file，因为不应该存在数据
    temp_log_file = tmp_path / "multi_block_test.log"
    # 创建超长记录（需要跨块存储）
    long_data = "A" * (LEVELDBLOG_BLOCK_LEN * 2)

    ds = DataStore()
    ds.open_for_write(str(temp_log_file))
    ds.write(long_data)
    ds.close()

    # 篡改第二个块中的数据
    with open(temp_log_file, "r+b") as f:
        # 定位到第二个块的数据部分
        f.seek(LEVELDBLOG_BLOCK_LEN + LEVELDBLOG_HEADER_LEN)
        f.write(b"TAMPERED")  # 篡改数据

    # 验证读取时触发异常
    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))
    with pytest.raises(ValidationError, match="Invalid record checksum"):
        next(ds)


def test_corrupted_header(temp_log_file):
    """损坏的文件头应触发异常"""
    with open(temp_log_file, "r+b") as f:
        f.seek(0)
        f.write(b"XXXX")  # 破坏标识符

    ds = DataStore()
    with pytest.raises(Exception, match="Invalid header"):
        ds.open_for_scan(str(temp_log_file))


def test_valid_records_after_tampered_one(temp_log_file):
    """验证篡改记录后的记录仍可正常读取"""
    # 篡改第一条记录
    with open(temp_log_file, "r+b") as f:
        f.seek(LEVELDBLOG_HEADER_LEN + 7)
        f.write(b"Tampered!!!!")

    ds = DataStore()
    ds.open_for_scan(str(temp_log_file))

    # 第一条记录应失败
    with pytest.raises(ValidationError):
        next(ds)

    # 第二条记录应正常读取
    record = next(ds)
    assert record == "Another valid record"

    # 后续无记录
    with pytest.raises(StopIteration):
        next(ds)
