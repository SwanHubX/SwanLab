"""
@author: cunyue
@file: test_system_hardware_e2e.py
@time: 2026/3/31
@description: 硬件采集流水线的测试，测试 new(ctx) 的完整编排逻辑
"""

from unittest.mock import patch

import pytest

from swanlab.sdk.internal.context import RunConfig, RunContext
from swanlab.sdk.internal.run.system import _new_raw
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings.metadata import EnvironmentSettings, MonitorSettings
from swanlab.sdk.typings.run.system import (
    AppleSiliconSnapshot,
    CPUSnapshot,
    GitSnapshot,
    MemorySnapshot,
    SystemEnvironment,
)

# ──────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────


@pytest.fixture
def make_ctx(tmp_path):
    """构建测试用 RunContext"""

    def _make(hardware=True, monitor_enable=True, git=True):
        settings = Settings(
            mode="disabled",
            environment=EnvironmentSettings(
                hardware=hardware,
                runtime=False,
                requirements=False,
                conda=False,
                git=git,
            ),
            monitor=MonitorSettings(enable=monitor_enable),
        )
        config = RunConfig(run_dir=tmp_path, settings=settings)
        return RunContext(config)

    return _make


def _mock_env_none():
    """返回隔离环境采集的 patch 上下文列表"""
    return [
        patch("swanlab.sdk.internal.run.system.git.get", return_value=None),
        patch("swanlab.sdk.internal.run.system.runtime.get", return_value=None),
        patch("swanlab.sdk.internal.run.system.conda.get", return_value=None),
        patch("swanlab.sdk.internal.run.system.requirements.get", return_value=None),
    ]


# ──────────────────────────────────────────────
# 1. Happy Path: 各平台正常采集
# ──────────────────────────────────────────────


class TestHardwareCollectionHappyPath:
    """测试不同平台下硬件信息的正常采集"""

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "darwin")
    def test_apple_silicon_macos(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """macOS Apple Silicon: apple_silicon 有值，cpu/memory 为 None"""
        mock_apple.get.return_value = AppleSiliconSnapshot(
            name="Apple M3 Pro", memory=18, memory_unit="GB", cpu_count=12
        )
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.apple_silicon is not None
            assert env.metadata.hardware.apple_silicon.name == "Apple M3 Pro"
            assert env.metadata.hardware.cpu is None
            assert env.metadata.hardware.memory is None
            assert env.shim.slug == "macos-arm"
            assert env.shim.enable_cpu is True
            assert env.shim.enable_memory is True
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_linux_cpu_memory(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """Linux: cpu 和 memory 有值，apple_silicon 为 None"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="Intel Core i7", physical_count=8, logical_count=16)
        mock_memory.get.return_value = MemorySnapshot(total=32, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.cpu is not None
            assert env.metadata.hardware.cpu.brand == "Intel Core i7"
            assert env.metadata.hardware.memory is not None
            assert env.metadata.hardware.memory.total == 32
            assert env.metadata.hardware.apple_silicon is None
            assert env.shim.slug == "linux"
            assert env.shim.enable_cpu is True
            assert env.shim.enable_memory is True
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "darwin")
    def test_macos_intel_fallback(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """macOS Intel: Apple.get() 返回 None，回退到 CPU + Memory"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="Intel Core i7", physical_count=4, logical_count=8)
        mock_memory.get.return_value = MemorySnapshot(total=16, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.apple_silicon is None
            assert env.metadata.hardware.cpu is not None
            assert env.metadata.hardware.memory is not None
            assert env.shim.slug == "macos-intel"
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "win32")
    def test_windows_cpu_memory(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """Windows: cpu 和 memory 有值"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="AMD Ryzen 9", physical_count=16, logical_count=32)
        mock_memory.get.return_value = MemorySnapshot(total=64, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.shim.slug == "windows"
            assert env.metadata.hardware.cpu is not None
            assert env.metadata.hardware.memory is not None
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "freebsd13")
    def test_unknown_platform(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """未知平台: 所有 vendor 返回 None，非 linux/darwin 平台 slug 为 windows"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = None
        mock_memory.get.return_value = None
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.shim.slug == "windows"  # 非 linux/darwin 均归为 windows
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.apple_silicon is None
            assert env.metadata.hardware.cpu is None
            assert env.metadata.hardware.memory is None
            assert env.shim.enable_cpu is False
            assert env.shim.enable_memory is False
        finally:
            for p in env_patches:
                p.stop()


# ──────────────────────────────────────────────
# 2. Settings 控制硬件采集
# ──────────────────────────────────────────────


class TestHardwareCollectionSettingsGating:
    """测试 Settings 对硬件采集行为的控制"""

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_both_disabled(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """hardware=False, monitor=False: metadata.hardware 为 None"""
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx(hardware=False, monitor_enable=False))
            assert env.metadata.hardware is None
        finally:
            for p in env_patches:
                p.stop()
        mock_apple.get.assert_not_called()
        mock_cpu.get.assert_not_called()
        mock_memory.get.assert_not_called()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_hardware_false_monitor_true(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """hardware=False, monitor=True: 内部采集了硬件（用于 shim），但 metadata 中被删除"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="Test CPU", physical_count=4, logical_count=8)
        mock_memory.get.return_value = MemorySnapshot(total=16, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx(hardware=False, monitor_enable=True))
            # del_hardware 被调用，metadata 中不包含硬件信息
            assert env.metadata.hardware is None
            # 但 shim 仍然基于硬件信息构建（用于监控）
            assert env.shim.slug == "linux"
            assert env.shim.enable_cpu is True
            assert env.shim.enable_memory is True
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_hardware_true_monitor_false(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """hardware=True, monitor=False: 正常采集硬件信息"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="Test CPU", physical_count=4, logical_count=8)
        mock_memory.get.return_value = MemorySnapshot(total=16, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx(hardware=True, monitor_enable=False))
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.cpu is not None
        finally:
            for p in env_patches:
                p.stop()


# ──────────────────────────────────────────────
# 3. 错误处理: vendor 采集失败或部分失败
# ──────────────────────────────────────────────


class TestHardwareCollectionErrorHandling:
    """测试 vendor 采集失败时的优雅降级"""

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_all_vendors_none(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """所有 vendor 返回 None: HardwareSnapshot 仍创建，所有字段为 None"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = None
        mock_memory.get.return_value = None
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.apple_silicon is None
            assert env.metadata.hardware.cpu is None
            assert env.metadata.hardware.memory is None
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_cpu_none_memory_ok(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """Apple=None, CPU=None, Memory 有效: cpu 为 None，memory 有值"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = None
        mock_memory.get.return_value = MemorySnapshot(total=32, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.cpu is None
            assert env.metadata.hardware.memory is not None
            assert env.metadata.hardware.memory.total == 32
            assert env.shim.enable_memory is True
            assert env.shim.enable_cpu is False
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_cpu_ok_memory_none(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """Apple=None, CPU 有效, Memory=None: cpu 有值，memory 为 None"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="AMD Ryzen 5", physical_count=6, logical_count=12)
        mock_memory.get.return_value = None
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.cpu is not None
            assert env.metadata.hardware.memory is None
            assert env.shim.enable_cpu is True
            assert env.shim.enable_memory is False
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_partial_environment(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """部分环境信息采集: git 有效，其余 None"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = None
        mock_memory.get.return_value = None
        git_snapshot = GitSnapshot(remote_url="https://github.com/test/repo", branch="main", commit="abc123")
        with (
            patch("swanlab.sdk.internal.run.system.git.get", return_value=git_snapshot),
            patch("swanlab.sdk.internal.run.system.runtime.get", return_value=None),
            patch("swanlab.sdk.internal.run.system.conda.get", return_value=None),
            patch("swanlab.sdk.internal.run.system.requirements.get", return_value=None),
        ):
            env, monitor = _new_raw(make_ctx(hardware=True, git=True))
            assert env.metadata.git is not None
            assert env.metadata.git.branch == "main"
            assert env.metadata.runtime is None
            assert env.conda is None
            assert env.requirements is None


# ──────────────────────────────────────────────
# 4. 数据一致性
# ──────────────────────────────────────────────


class TestHardwareSnapshotConsistency:
    """测试返回数据的结构一致性"""

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_return_type(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """new(ctx) 返回 (SystemEnvironment, None)"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = None
        mock_memory.get.return_value = None
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            result = _new_raw(make_ctx())
            assert isinstance(result, tuple)
            assert len(result) == 2
            assert isinstance(result[0], SystemEnvironment)
            assert result[1] is not None
        finally:
            for p in env_patches:
                p.stop()

    @patch("swanlab.sdk.internal.run.system.Memory")
    @patch("swanlab.sdk.internal.run.system.CPU")
    @patch("swanlab.sdk.internal.run.system.Apple")
    @patch("swanlab.sdk.internal.run.system.sys.platform", "linux")
    def test_accelerators_always_empty(self, mock_apple, mock_cpu, mock_memory, make_ctx):
        """accelerators 始终为空列表（尚未实现）"""
        mock_apple.get.return_value = None
        mock_cpu.get.return_value = CPUSnapshot(brand="Test", physical_count=4, logical_count=8)
        mock_memory.get.return_value = MemorySnapshot(total=16, total_unit="GB")
        env_patches = _mock_env_none()
        for p in env_patches:
            p.start()
        try:
            env, monitor = _new_raw(make_ctx())
            assert env.metadata.hardware is not None
            assert env.metadata.hardware.accelerators == []
            assert env.shim.accelerators == []
        finally:
            for p in env_patches:
                p.stop()
