"""
@author: cunyue
@file: test_log.py
@time: 2026/3/10
@description: 测试 swanlab.sdk.pkg.log 诊断日志模块 (极简黑盒版)
"""

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
    """测试暴露的日志格式符合预期"""
    log_dir = tmp_path / "logs"
    log_dir.mkdir()
    log_mod.bindfile(log_dir)

    log_mod.error("check format")

    content = (log_dir / "debug.log").read_text(encoding="utf-8")
    # 预期格式：时间 | ERROR   | swanlab.internal | 消息
    assert " | ERROR   | swanlab.internal | check format" in content
