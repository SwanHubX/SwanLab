"""
@author: caddiesnew
@file: base.py
@time: 2026/4/20
@description: 所有实体类的公共基类
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterator, Optional

from swanlab.sdk.internal.pkg import safe

from .typings.common import ApiResponseType

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
    def to_dict(self) -> Dict[str, Any]:
        """将实体序列化为 JSON 可序列化的字典。"""

    def _safe_request(self, method: Callable, path: str, **kwargs) -> ApiResponseType:
        """安全请求包装：捕获所有异常，始终返回 ApiResponse 而不抛出。"""

        def _on_error(e: BaseException) -> None:
            _err_msg[0] = str(e)

        _err_msg: list[Optional[str]] = [None]
        with safe.block(message=f"API request failed: {path}", on_error=_on_error):
            data = method(path, **kwargs).data
            return ApiResponseType(ok=True, data=data)
        result = ApiResponseType(ok=False, errmsg=_err_msg[0] or "request failed")
        self._errors.append(result.errmsg)
        return result

    def _get(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.get, path, **kwargs)

    def _post(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.post, path, **kwargs)

    def _put(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.put, path, **kwargs)

    def _delete(self, path: str, **kwargs) -> ApiResponseType:
        return self._safe_request(self._client.delete, path, **kwargs)

    def _build_url(self, path: str) -> str:
        return f"{self._api_host}/{path}"

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
            body = resp.data
            items = body.get("list", []) if isinstance(body, dict) else body
            if not items:
                break
            yield from items
            total_pages = body.get("pages", 1) if isinstance(body, dict) else 1
            if page >= total_pages:
                break
            page += 1

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        ident = getattr(self, "_path", None) or getattr(self, "_username", None) or "?"
        return f"{cls}('{ident}')"
