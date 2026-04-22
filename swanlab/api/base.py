"""
@author: caddiesnew
@file: base.py
@time: 2026/4/20
@description: 所有实体类的公共基类
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterator, Optional

from swanlab.sdk.internal.pkg import safe

from .typings.common import ApiPaginationType, ApiResponseType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class BaseEntity(ABC):
    """
    swanlab/api 实体类公共基类。

    统一持有 _client、_web_host 和 _api_host，提供 _get/_post/_put/_delete HTTP 快捷方法和 _paginate 分页迭代。
    所有 HTTP 请求通过 _safe_request 包裹，保证任何异常都不会导致程序 crash，统一返回 ApiResponse。
    子类只需实现 to_dict() 和业务逻辑。
    """

    def __init__(self, client: "Client", web_host: str, api_host: str) -> None:
        self._client: "Client" = client
        self._web_host: str = web_host
        self._api_host: str = api_host
        self._errors: list[str] = []

    @abstractmethod
    def json(self) -> Dict[str, Any]:
        """将实体序列化为 JSON 可序列化的字典。"""

    def _safe_request(self, method: Callable, path: str, **kwargs) -> ApiResponseType:
        """安全请求包装：捕获所有异常，始终返回 ApiResponse 而不抛出。"""
        _err: list[str] = []
        common_err: str = f"API request failed: {path}"

        @safe.decorator(message=common_err, on_error=lambda e: _err.append(str(e)))
        def _do():
            return method(path, **kwargs).data

        data = _do()
        if data is not None:
            return ApiResponseType(ok=True, data=data)
        errmsg = _err[0] if _err else common_err
        self._errors.append(errmsg)
        return ApiResponseType(ok=False, errmsg=errmsg)

    def _get(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.get, path, **kwargs)

    def _post(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.post, path, **kwargs)

    def _put(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.put, path, **kwargs)

    def _delete(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.delete, path, **kwargs)

    def _build_web_url(self, path: str) -> str:
        """构建前端 Web 页面 URL（使用 _web_host 而非 _api_host）。"""
        return f"{self._web_host}/{path}"

    def _paginate(self, path: str, *, page_size: int = 20, params: Optional[dict] = None) -> Iterator[dict]:
        """通用分页迭代器，自动处理 page/size 参数。"""
        page = 1
        while True:
            p = {"page": page, "size": page_size}
            if params:
                p.update({k: v for k, v in params.items() if v is not None})
            resp = self._get(path, params=p)
            if not resp.ok:
                return
            body: ApiPaginationType = resp.data
            items = body.get("list", [])
            if not items:
                break
            yield from items
            if page >= body.get("pages", 1):
                break
            page += 1

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        ident = getattr(self, "_path", None) or getattr(self, "_username", None) or "?"
        return f"{cls}('{ident}')"
