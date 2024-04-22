#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-22 18:23:47
@File: swanlab\intergration\utils\autologging.py
@IDE: vscode
@Description:
    autologging,用于实现集成指标的获取与可视化
"""
import asyncio
import functools
import inspect
import logging
import sys
from typing import Any, Dict, Optional, Sequence, TypeVar


if sys.version_info >= (3, 8):
    from typing import Protocol
else:
    from typing_extensions import Protocol


K = TypeVar("K", bound=str)
V = TypeVar("V")


class Response(Protocol[K, V]):
    def __getitem__(self, key: K) -> V: ...  # pragma: no cover

    def get(self, key: K, default: Optional[V] = None) -> Optional[V]: ...  # pragma: no cover


class ArgumentResponseResolver(Protocol):
    def __call__(
        self,
        args: Sequence[Any],
        kwargs: Dict[str, Any],
        response: Response,
        start_time: float,
        time_elapsed: float,
    ) -> Optional[Dict[str, Any]]: ...  # pragma: no cover
