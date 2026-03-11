"""
@author: cunyue
@file: test_fs_write.py
@time: 2026/3/11 14:48
@description: 测试 SwanLab SDK 文件系统辅助函数写入功能
"""

import io
from pathlib import Path
from unittest.mock import patch

import pytest

# 假设你的模块名为 write，根据实际路径调整 import
from swanlab.sdk.internal.pkg.fs.write import safe_write


def test_safe_write_to_io_text(tmp_path: Path):
    """测试写入到已打开的文件句柄 (文本模式)"""
    file_path = tmp_path / "test_io.txt"

    with open(file_path, "w", encoding="utf-8") as f:
        safe_write(f, "hello swanlab")

    assert file_path.read_text(encoding="utf-8") == "hello swanlab"


def test_safe_write_to_io_bytes(tmp_path: Path):
    """测试写入到已打开的 BytesIO 对象 (二进制模式，不带 fileno)"""
    # 使用 BytesIO 测试，因为它没有 fileno 属性，可以测试 getattr 容错逻辑
    mem_file = io.BytesIO()
    safe_write(mem_file, b"\x00\x01\x02")

    assert mem_file.getvalue() == b"\x00\x01\x02"


def test_safe_write_path_non_atomic_text_default_encoding(tmp_path: Path):
    """测试非原子写入：普通文本，默认 utf-8 编码"""
    file_path = tmp_path / "test_text_default.txt"
    safe_write(file_path, "默认编码测试")

    assert file_path.exists()
    assert file_path.read_text(encoding="utf-8") == "默认编码测试"


def test_safe_write_path_non_atomic_text_custom_encoding(tmp_path: Path):
    """测试非原子写入：指定自定义编码 (如 GBK)"""
    file_path = tmp_path / "test_text_gbk.txt"
    # 使用 gbk 写入
    safe_write(file_path, "GBK编码测试", encoding="gbk")

    assert file_path.exists()
    # 验证确实是用 GBK 写入的，如果用 utf-8 读会报错或乱码
    assert file_path.read_text(encoding="gbk") == "GBK编码测试"


def test_safe_write_path_non_atomic_binary(tmp_path: Path):
    """测试非原子写入：二进制数据"""
    file_path = tmp_path / "test_bin.bin"
    # 传入二进制内容和 wb 模式
    safe_write(file_path, b"binary content", mode="wb")

    assert file_path.exists()
    assert file_path.read_bytes() == b"binary content"


def test_safe_write_path_atomic_text(tmp_path: Path):
    """测试原子写入：多级目录，文本模式"""
    # 故意使用多级不存在的目录，顺便测试 safe_mkdir 是否被正确调用
    file_path = tmp_path / "nested" / "dir" / "test_atomic.txt"
    safe_write(file_path, "atomic content", atomic=True)

    assert file_path.read_text(encoding="utf-8") == "atomic content"
    # 验证同级目录下没有残留的 .tmp_swanlab_ 隐藏文件
    parent_dir = file_path.parent
    assert len(list(parent_dir.iterdir())) == 1
    assert list(parent_dir.iterdir())[0].name == "test_atomic.txt"


def test_safe_write_path_atomic_binary(tmp_path: Path):
    """测试原子写入：二进制模式"""
    file_path = tmp_path / "test_atomic.bin"
    safe_write(file_path, b"\x01\x02\x03", mode="wb", atomic=True)

    assert file_path.read_bytes() == b"\x01\x02\x03"


@patch("swanlab.sdk.internal.pkg.fs.write.os.replace")
def test_safe_write_atomic_rollback(mock_replace, tmp_path: Path):
    """
    【灾难恢复测试】测试原子写入在 os.replace 阶段遭遇异常（如权限问题、被其他进程锁死），
    是否能够成功执行清理逻辑，不留垃圾。
    """
    file_path = tmp_path / "test_rollback.txt"

    # 模拟 os.replace 时突然发生异常
    mock_replace.side_effect = PermissionError("Simulated Replace Error")

    with pytest.raises(PermissionError, match="Simulated Replace Error"):
        safe_write(file_path, "some content that will fail", atomic=True)

    # 核心断言：异常发生后，临时文件被 os.remove(temp_path) 清理掉了
    # 此时 tmp_path 应该是一个空目录，不留任何 .tmp_swanlab_ 文件
    assert len(list(tmp_path.iterdir())) == 0
    assert not file_path.exists()
