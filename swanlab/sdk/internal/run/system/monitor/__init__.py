"""
@author: cunyue
@file: __init__.py
@time: 2026/4/9 16:41
@description: SwanLab 硬件监控对象
"""

from typing import TYPE_CHECKING, Optional

from swanlab.sdk.internal.pkg.timer import Timer
from swanlab.sdk.typings.run.system import SystemShim

if TYPE_CHECKING:
    from swanlab import Run

__all__ = ["Monitor"]


class Monitor:
    def __init__(self, shim: SystemShim):
        self._timer: Optional[Timer] = None
        self._shim = shim

    @property
    def is_running(self):
        return self._timer is not None and self._timer.is_running

    def start(self, run: "Run"):
        pass

    def cancel(self) -> None:
        assert self._timer is not None, "HardwareMonitor is not running"
        self._timer.cancel()

    def join(self) -> None:
        assert self._timer is not None, "HardwareMonitor is not running"
        self._timer.join()
