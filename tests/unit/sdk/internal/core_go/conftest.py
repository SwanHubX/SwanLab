"""core_go 测试共享夹具与 Fake 定义。"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import pytest

REPO_ROOT = Path(__file__).resolve().parents[5]


class FakeLauncher:
    """CoreLauncher 共享测试替身，支持 owner/reuse/spinup 场景配置与调用记录。"""

    def __init__(
        self,
        binary: Optional[Path] = None,
        *,
        detach: bool = False,
        idle_timeout: float = 0.0,
        owns_process: bool = True,
        exit_code: Optional[int] = None,
        channel: Any = None,
        start_error: Optional[BaseException] = None,
    ) -> None:
        self.binary = binary
        self.detach = detach
        self.idle_timeout = idle_timeout
        self.owns_process = owns_process
        self.process_exited = False
        self.channel = object() if channel is None else channel
        self.started = False
        self.starts = 0
        self.closed = 0
        self.aborted = 0
        self.exit_code = exit_code
        self.start_error = start_error

    def start(self) -> None:
        self.started = True
        self.starts += 1
        if self.start_error is not None:
            raise self.start_error

    def close(self) -> None:
        self.closed += 1

    def abort(self) -> None:
        self.aborted += 1


@pytest.fixture
def short_tmp_path(tmp_path):
    """返回足够短的临时目录，避免 Unix socket 路径超限。

    POSIX 上 Unix domain socket 路径通常限制在 104 字节以内，pytest 的
    ``tmp_path`` 可能超过该限制；Windows 仍使用 ``tmp_path``。
    """
    if sys.platform == "win32":
        yield tmp_path
        return

    root = Path(tempfile.mkdtemp(prefix="swl-t-", dir="/tmp"))
    try:
        yield root
    finally:
        shutil.rmtree(root, ignore_errors=True)


@pytest.fixture
def fake_launcher_cls():
    """提供可直连构造的 FakeLauncher 类。"""
    return FakeLauncher


@pytest.fixture
def install_fake_launcher(monkeypatch):
    """在目标模块安装并记录 FakeLauncher 实例与启动序号的夹具工厂。"""

    def _install(
        target_module: Any,
        *,
        owns_process: Union[bool, Callable[[int], bool]] = True,
        exit_code: Union[Optional[int], Callable[[bool], Optional[int]]] = None,
        channel_factory: Optional[Callable[[int, bool], Any]] = None,
        abort_increments_closed: bool = False,
    ) -> Tuple[List[FakeLauncher], List[int]]:
        instances: List[FakeLauncher] = []
        starts: List[int] = []

        class RecordingLauncher(FakeLauncher):
            def __init__(
                self, binary: Optional[Path] = None, *, detach: bool = False, idle_timeout: float = 0.0
            ) -> None:
                index = len(instances)
                owned = owns_process(index) if callable(owns_process) else owns_process
                code = exit_code(detach) if callable(exit_code) else exit_code
                channel = channel_factory(index, detach) if channel_factory is not None else f"ch-{index}"
                super().__init__(
                    binary,
                    detach=detach,
                    idle_timeout=idle_timeout,
                    owns_process=owned,
                    exit_code=code,
                    channel=channel,
                )
                self.index = index
                instances.append(self)

            def start(self) -> None:
                starts.append(self.index)
                super().start()

            def abort(self) -> None:
                super().abort()
                if abort_increments_closed:
                    self.closed += 1

        monkeypatch.setattr(target_module, "CoreLauncher", RecordingLauncher)
        return instances, starts

    return _install


@pytest.fixture
def run_python_isolated():
    """以仓库根路径在干净子进程中执行一段 Python 导入/探测脚本。"""

    def _run(
        code: str,
        *,
        env_overrides: Optional[Dict[str, str]] = None,
        drop_env: Sequence[str] = (),
        check: bool = True,
    ) -> subprocess.CompletedProcess:
        env = {k: v for k, v in os.environ.items() if k not in drop_env}
        env["PYTHONPATH"] = str(REPO_ROOT)
        if env_overrides:
            env.update(env_overrides)
        return subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(REPO_ROOT),
            env=env,
            capture_output=True,
            text=True,
            check=check,
        )

    return _run
