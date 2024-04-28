#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-23 11:52:13
@File: swanlab\utils\timer.py
@IDE: vscode
@Description:
    计时器
"""

import time
from typing import Any


class Timer:
    def __init__(self) -> None:
        self.start_time: float = time.time()
        self.start: float = time.perf_counter()
        self.stop: float = self.start

    def __enter__(self) -> "Timer":
        return self

    def __exit__(self, *args: Any) -> None:
        self.stop = time.perf_counter()

    @property
    def elapsed(self) -> float:
        return self.stop - self.start
