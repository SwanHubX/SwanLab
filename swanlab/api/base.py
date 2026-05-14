"""
@author: caddiesnew
@file: base.py
@time: 2026/4/20
@description: 所有实体类的公共基类
"""

import asyncio
import random
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple

from swanlab.sdk.internal.pkg import safe

from .typings.common import MAX_CONCURRENT_COUNT, ApiPaginationType, ApiResponseType, PaginatedQuery

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


@dataclass(frozen=True)
class ApiClientContext:
    """Shared context: all sub-entities reuse the same authenticated client instance."""

    client: "Client"
    web_host: str
    api_host: str
    username: str
    name: str


class BaseEntity(ABC):
    """
    Base class for all swanlab/api entities.

    Holds a shared ``_ctx`` (ApiClientContext) and provides ``_get``/``_post``/``_put``/``_delete``
    HTTP shortcuts as well as ``_paginate`` for paginated iteration.
    All HTTP requests are wrapped by ``_safe_request``, which guarantees that no exception
    crashes the process — every call returns an ``ApiResponseType``.
    Subclasses only need to implement ``json()`` and their business logic.
    """

    def __init__(self, ctx: ApiClientContext) -> None:
        self._ctx: ApiClientContext = ctx
        self._errors: list[str] = []

    def _ensure_data(self) -> Any:
        """Lazy-load data. Override in single-entity subclasses; iterator subclasses need not."""
        return None

    def wrapper(self) -> ApiResponseType:
        """Eager mode: trigger ``_ensure_data`` and return ``ApiResponseType`` based on ``_errors``."""
        self._ensure_data()
        if self._errors:
            return ApiResponseType(ok=False, errmsg=self._errors[-1])
        return ApiResponseType(ok=True, data=self)

    @abstractmethod
    def json(self) -> Dict[str, Any]:
        """Serialize the entity to a JSON-compatible dict."""

    def _safe_request(self, method: Callable, path: str, **kwargs) -> ApiResponseType:
        """Safe request wrapper: catches all exceptions and always returns ``ApiResponseType``."""
        _err: list[str] = []
        common_err: str = f"API Request Failed: {path}"

        @safe.decorator(message=None, on_error=lambda e: _err.append(str(e)))
        def _do():
            return method(path, **kwargs).data

        data = _do()
        if data is not None:
            return ApiResponseType(ok=True, data=data)
        errmsg = _err[0] if _err else common_err
        self._errors.append(errmsg)
        return ApiResponseType(ok=False, errmsg=errmsg)

    def _get(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._ctx.client.get, path, **kwargs)

    def _post(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._ctx.client.post, path, **kwargs)

    def _put(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._ctx.client.put, path, **kwargs)

    def _delete(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._ctx.client.delete, path, **kwargs)

    def _concurrent_request(
        self,
        requests: Sequence[Tuple[Callable, str, Dict[str, Any]]],
        semaphore: int = MAX_CONCURRENT_COUNT,
    ) -> List[ApiResponseType]:
        """
        Execute multiple HTTP requests concurrently, bounded by a semaphore.

        :param requests: List of ``(method, path, kwargs)`` tuples.
            ``method``: one of ``self._get`` / ``self._post`` / ``self._put`` / ``self._delete``.
            ``path``: API endpoint path.
            ``kwargs``: request parameters (``params=``, ``data=``, etc.).
        :param semaphore: Maximum number of concurrent requests.
        :return: ``ApiResponseType`` list with the same length as *requests*.
        """
        if not requests:
            return []

        async def _run():
            sem = asyncio.Semaphore(semaphore)

            async def _fetch(req: Tuple[Callable, str, Dict[str, Any]]) -> ApiResponseType:
                method, path, kwargs = req
                async with sem:
                    loop = asyncio.get_running_loop()
                    return await loop.run_in_executor(None, lambda: method(path, **kwargs))

            return await asyncio.gather(*[_fetch(r) for r in requests])

        return asyncio.run(_run())

    def _build_web_url(self, path: str) -> str:
        """Build a frontend web page URL (uses ``_web_host`` instead of ``_api_host``)."""
        return f"{self._ctx.web_host}/{path}"

    def _paginate(
        self,
        path: str,
        query: PaginatedQuery,
        *,
        page_info: Dict[str, Any],
        extra: Optional[Dict[str, Any]] = None,
    ) -> Iterator[dict]:
        """Generic paginated iterator driven by ``PaginatedQuery``."""
        page = query.page
        while True:
            p = query.to_params(**(extra or {}))
            # 覆盖当前页码（翻页时自增）
            p["page"] = page
            resp = self._get(path, params=p)
            if not resp.ok:
                return
            body: ApiPaginationType = resp.data

            if page_info["total"] == 0:
                page_info["total"] = body.get("total", 0)
            if page_info["pages"] == 0:
                page_info["pages"] = body.get("pages", 1)
            items = body.get("list", [])
            if not items:
                break
            yield from items
            # 随机休眠控制 qps
            time.sleep(random.random())
            if page >= body.get("pages", 1):
                break
            if not query.all:
                break
            page += 1

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        ident = getattr(self, "_path", None) or getattr(self, "_username", None) or "?"
        return f"{cls}('{ident}')"
