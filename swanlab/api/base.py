"""
@author: caddiesnew
@file: base.py
@time: 2026/4/20
@description: 所有实体类的公共基类
"""

import random
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterator, Optional

from swanlab.sdk.internal.pkg import safe

from .typings.common import ApiPaginationType, ApiResponseType, PaginatedQuery

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


@dataclass(frozen=True)
class ApiClientContext:
    """共享上下文：所有子实体复用同一个登录态实例。"""

    client: "Client"
    web_host: str
    api_host: str
    username: str
    name: str


class BaseEntity(ABC):
    """
    swanlab/api 实体类公共基类。

    统一持有 _ctx（_ApiContext），提供 _get/_post/_put/_delete HTTP 快捷方法和 _paginate 分页迭代。
    所有 HTTP 请求通过 _safe_request 包裹，保证任何异常都不会导致程序 crash，统一返回 ApiResponse。
    子类只需实现 json() 和业务逻辑。
    """

    def __init__(self, ctx: ApiClientContext) -> None:
        self._ctx: ApiClientContext = ctx
        self._errors: list[str] = []

    def _ensure_data(self) -> Any:
        """按需加载数据。单实体子类重写此方法；迭代器子类无需重写。"""
        return None

    def wrapper(self) -> ApiResponseType:
        """Eager 模式：触发子类 _ensure_data 加载数据，根据 _errors 返回 ApiResponseType。"""
        self._ensure_data()
        if self._errors:
            return ApiResponseType(ok=False, errmsg=self._errors[-1])
        return ApiResponseType(ok=True, data=self)

    @abstractmethod
    def json(self) -> Dict[str, Any]:
        """将实体序列化为 JSON 可序列化的字典。"""

    def _safe_request(self, method: Callable, path: str, **kwargs) -> ApiResponseType:
        """安全请求包装：捕获所有异常，始终返回 ApiResponseType 而不抛出。"""
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

    def _build_web_url(self, path: str) -> str:
        """构建前端 Web 页面 URL（使用 _web_host 而非 _api_host）。"""
        return f"{self._ctx.web_host}/{path}"

    def _paginate(
        self,
        path: str,
        query: PaginatedQuery,
        *,
        page_info: Dict[str, Any],
        extra: Optional[Dict[str, Any]] = None,
    ) -> Iterator[dict]:
        """通用分页迭代器，基于 PaginatedQuery 驱动翻页逻辑。"""
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
