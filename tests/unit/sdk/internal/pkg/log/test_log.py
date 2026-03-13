"""
@author: cunyue
@file: test_log.py
@time: 2026/3/10
@description: 测试 swanlab.sdk.pkg.log 诊断日志模块 (极简黑盒版)
"""

import os
import stat

import pytest

import swanlab.sdk.internal.pkg.log as log_mod


def test_log_lifecycle_and_output(tmp_path):
    """测试核心业务流：缓冲阶段记录 -> 绑定文件触发落盘 -> 直接写入文件"""
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    # 1. 缓冲阶段（尚未 bindfile）
    log_mod.debug("buffered debug")
    log_mod.info("buffered info")

    # 2. 绑定文件（之前的缓冲应自动落盘）
    log_mod.bindfile(log_dir)

    # 3. 绑定后写入（应直接落盘）
    log_mod.warning("direct warning")
    log_mod.error("direct error")
    log_mod.critical("direct critical")

    # 验证最终落盘内容
    log_file = log_dir / "debug.log"
    content = log_file.read_text(encoding="utf-8")

    assert "buffered debug" in content
    assert "buffered info" in content
    assert "Diagnostic log bound to file" in content
    assert "direct warning" in content
    assert "direct error" in content
    assert "direct critical" in content


def test_bindfile_invalid_dir(tmp_path):
    """测试绑定到不存在的目录时抛出 FileNotFoundError"""
    with pytest.raises(FileNotFoundError, match="Log directory does not exist"):
        log_mod.bindfile(tmp_path / "not_exist")


def test_bindfile_is_idempotent(tmp_path):
    """测试重复调用 bindfile 的幂等性（不报错且不重复创建资源）"""
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    # 第一次绑定
    log_mod.bindfile(log_dir)
    log_mod.info("first bound msg")

    # 第二次重复绑定（应静默忽略）
    log_mod.bindfile(log_dir)
    log_mod.info("second bound msg")

    content = (log_dir / "debug.log").read_text(encoding="utf-8")

    assert "first bound msg" in content
    assert "second bound msg" in content


def test_log_format(tmp_path):
    """测试写入日志文件的内容与传入消息一致（formatter 为纯消息模式）"""
    log_dir = tmp_path / "logs"
    log_dir.mkdir()
    log_mod.bindfile(log_dir)

    log_mod.error("check format")

    content = (log_dir / "debug.log").read_text(encoding="utf-8")
    assert "check format" in content


def test_secure_file_permissions(tmp_path):
    """测试生成的诊断日志文件是否具备安全的 0600 权限 (仅限 POSIX 系统)"""
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    # 绑定并写入一条日志，确保底层 _open 方法被真实触发
    log_mod.bindfile(log_dir)
    log_mod.info("trigger file creation for permission check")

    log_file = log_dir / "debug.log"
    assert log_file.exists()

    # 仅在类 Unix 系统下校验 POSIX 权限
    if os.name == "posix":
        # 获取文件的完整状态信息
        mode = os.stat(log_file).st_mode
        # 提取底部的八进制权限位
        permissions = stat.S_IMODE(mode)

        # 验证权限是否被严格限制为 0o600 (即 rw-------，仅所有者可读写)
        assert permissions == 0o600
