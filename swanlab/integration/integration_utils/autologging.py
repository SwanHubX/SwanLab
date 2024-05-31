#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-04-22 18:23:47
@File: swanlab\intergration\utils\autologging.py
@IDE: vscode
@Description:
    autologging,用于实现集成指标的获取
"""
import asyncio
import functools
import inspect
import swanlab
import sys
import time
from swanlab.data.run import SwanLabRun
from typing import Any, Dict, Optional, Sequence, TypeVar, TYPE_CHECKING
from swanlab.integration.integration_utils.get_modules import get_module
from swanlab.integration.integration_utils.timer import Timer


if sys.version_info >= (3, 8):
    from typing import Protocol
else:
    from typing_extensions import Protocol

if TYPE_CHECKING:
    import openai


AutologInitArgs = Optional[Dict[str, Any]]

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


class PatchAPI:
    def __init__(
        self, name: str, symbols: Sequence[str], resolver: ArgumentResponseResolver, client, lib_version
    ) -> None:
        """Patches the API to log SwanLab Media or metrics."""
        # name of the LLM provider, e.g. "Cohere" or "OpenAI" or package name like "Transformers"
        self.name = name
        # api library name or client name, e.g.  "openai" or "openai.OpenAI()"
        self._api = client
        # api library version
        self.lib_version = lib_version
        # dictionary of original methods
        self.original_methods: Dict[str, Any] = {}
        # list of symbols to patch, e.g. ["Client.generate", "Edit.create"] or ["Pipeline.__call__"]
        self.symbols = symbols
        # resolver callable to convert args/response into a dictionary of SwanLab media objects or metrics
        self.resolver = resolver

    @property
    def set_api(self) -> Any:
        """Returns the API module."""
        lib_name = self.name.lower()
        if self._api is None:
            self._api = get_module(
                name=lib_name,
                required=f"To use the SwanLab {self.name} Autolog, "
                f"you need to have the `{lib_name}` python "
                f"package installed. Please install it with `pip install {lib_name}`.",
                lazy=False,
            )
        return self._api

    def patch(self, run: "SwanLabRun") -> None:
        """Patches the API to log media or metrics to SwanLab."""
        for symbol in self.symbols:
            # split on dots, e.g. "Client.generate" -> ["Client", "generate"]
            symbol_parts = symbol.split(".")

            # and get the attribute from the module
            original = functools.reduce(getattr, symbol_parts, self.set_api)

            def method_factory(original_method: Any):
                async def async_method(*args, **kwargs):
                    future = asyncio.Future()

                    async def callback(coro):
                        try:
                            result = await coro
                            loggable_dict = self.resolver(
                                args, kwargs, result, self.lib_version, timer.start_time, time_elapsed
                            )
                            if loggable_dict is not None:
                                swanlab.log(loggable_dict)
                            future.set_result(result)
                        except Exception as e:
                            print(e)

                    with Timer() as timer:
                        start_time = time.perf_counter()
                        coro = original_method(*args, **kwargs)
                        end_time = time.perf_counter()
                        time_elapsed = int(end_time - start_time + 0.5)
                        asyncio.ensure_future(callback(coro))

                    return await future

                def sync_method(*args, **kwargs):
                    with Timer() as timer:
                        start_time = time.perf_counter()
                        result = original_method(*args, **kwargs)
                        end_time = time.perf_counter()
                        time_elapsed = int(end_time - start_time + 0.5)
                        try:
                            loggable_dict = self.resolver(
                                args,
                                kwargs,
                                result,
                                self.lib_version,
                                timer.start_time,
                                time_elapsed,
                            )
                            if loggable_dict is not None:
                                swanlab.log(loggable_dict)
                        except Exception as e:
                            print(e)
                        return result

                if inspect.iscoroutinefunction(original_method):
                    return functools.wraps(original_method)(async_method)
                else:
                    return functools.wraps(original_method)(sync_method)

            # save original method
            self.original_methods[symbol] = original
            # monkey patch the method
            if len(symbol_parts) == 1:
                setattr(self.set_api, symbol_parts[0], method_factory(original))
            else:
                setattr(
                    functools.reduce(getattr, symbol_parts[:-1], self.set_api),
                    symbol_parts[-1],
                    method_factory(original),
                )

    def unpatch(self) -> None:
        """Unpatches the API."""
        for symbol, original in self.original_methods.items():
            # split on dots, e.g. "Client.generate" -> ["Client", "generate"]
            symbol_parts = symbol.split(".")
            # unpatch the method
            if len(symbol_parts) == 1:
                setattr(self.set_api, symbol_parts[0], original)
            else:
                setattr(
                    functools.reduce(getattr, symbol_parts[:-1], self.set_api),
                    symbol_parts[-1],
                    original,
                )


class AutologAPI:
    def __init__(
        self,
        name: str,
        symbols: Sequence[str],
        resolver: ArgumentResponseResolver,
        client=None,
        lib_version=None,
        mode: str = None,
    ) -> None:
        """Autolog API calls to SwanLab."""

        self._patch_api = PatchAPI(
            name=name, symbols=symbols, resolver=resolver, client=client, lib_version=lib_version
        )
        self._name = self._patch_api.name
        self._run: Optional[SwanLabRun] = None
        self.mode = mode
        self.client: openai.Client = client
        self.lib_version = lib_version

    @property
    def _is_enabled(self) -> bool:
        """Returns whether autologging is enabled."""
        return self._run is not None

    def __call__(self, init: AutologInitArgs = None) -> None:
        """Enable autologging."""
        self.enable(init=init)

    def _run_init(self, init: AutologInitArgs = None) -> None:
        """Handle SwanLab run initialization."""
        # - autolog(init: dict = {...}) calls swanlab.init(**{...})
        #   we only track if the run was created by autolog
        #    - todo: autolog(init: dict | run = run) would use the user-provided run
        # - autolog() calls swanlab.init()
        if init:
            self._run = swanlab.init(**init)

    def enable(self, init: AutologInitArgs = None) -> None:
        """Enable autologging.

        Args:
            init: Optional dictionary of arguments to pass to SwanLab.init().

        """
        if self._is_enabled:
            print(f"{self._name} autologging is already enabled, disabling and re-enabling.")
            self.disable()

        print(f"Enabling {self._name} autologging.")
        self._run_init(init=init)

        self._patch_api.patch(self._run)

    def disable(self) -> None:
        """Disable autologging."""
        if self._run is None:
            return

        print(f"Disabling {self._name} autologging.")

        self._run = None

        self._patch_api.unpatch()


if __name__ == "__main__":
    print("This is the main program")
