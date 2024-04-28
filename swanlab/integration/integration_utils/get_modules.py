#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-22 20:10:01
@File: swanlab\utils\get_modules.py
@IDE: vscode
@Description:
    获取库名称
"""

import importlib
import importlib.util
import sys
import threading
import types
from dataclasses import asdict, is_dataclass
from datetime import date, datetime, timedelta
from importlib import import_module
from sys import getsizeof
from types import ModuleType
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Generator,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    TextIO,
    Tuple,
    TypeVar,
    Union,
)

_not_importable = set()


class LazyModuleState:
    def __init__(self, module: types.ModuleType) -> None:
        self.module = module
        self.load_started = False
        self.lock = threading.RLock()

    def load(self) -> None:
        with self.lock:
            if self.load_started:
                return
            self.load_started = True
            assert self.module.__spec__ is not None
            assert self.module.__spec__.loader is not None
            self.module.__spec__.loader.exec_module(self.module)
            self.module.__class__ = types.ModuleType


class LazyModule(types.ModuleType):
    def __getattribute__(self, name: str) -> Any:
        state = object.__getattribute__(self, "__lazy_module_state__")
        state.load()
        return object.__getattribute__(self, name)

    def __setattr__(self, name: str, value: Any) -> None:
        state = object.__getattribute__(self, "__lazy_module_state__")
        state.load()
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        state = object.__getattribute__(self, "__lazy_module_state__")
        state.load()
        object.__delattr__(self, name)


def import_module_lazy(name: str) -> types.ModuleType:
    """Import a module lazily, only when it is used.

    Inspired by importlib.util.LazyLoader, but improved so that the module loading is
    thread-safe. Circular dependency between modules can lead to a deadlock if the two
    modules are loaded from different threads.

    """
    try:
        return sys.modules[name]
    except KeyError:
        spec = importlib.util.find_spec(name)
        if spec is None:
            raise ModuleNotFoundError
        module = importlib.util.module_from_spec(spec)
        module.__lazy_module_state__ = LazyModuleState(module)  # type: ignore
        module.__class__ = LazyModule
        sys.modules[name] = module
        return module


def get_module(
    name: str,
    required: Optional[Union[str, bool]] = None,
    lazy: bool = True,
) -> Any:
    """Return module or None. Absolute import is required.

    :param (str) name: Dot-separated module path. E.g., 'scipy.stats'.
    :param (str) required: A string to raise a ValueError if missing
    :param (bool) lazy: If True, return a lazy loader for the module.
    :return: (module|None) If import succeeds, the module will be returned.
    """
    if name not in _not_importable:
        try:
            if not lazy:
                return import_module(name)
            else:
                return import_module_lazy(name)
        except Exception:
            _not_importable.add(name)
            msg = f"Error importing optional module {name}"
            if required:
                print(msg)


def get_optional_module(name) -> Optional["importlib.ModuleInterface"]:  # type: ignore
    return get_module(name)
