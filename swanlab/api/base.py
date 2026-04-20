"""
@author: caddiesnew
@file: base.py
@time: 2026/4/20
@description: 所有实体类的公共基类
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class BaseEntity(ABC):
    """
    swanlab/api 实体类公共基类。

    统一持有 _client 和 _web_host，提供 _get/_post/_put/_delete HTTP 快捷方法和 _paginate 分页迭代。
    子类只需实现 to_dict() 和业务逻辑。
    """

    def __init__(self, client: "Client", web_host: str, api_host: str) -> None:
        self._client: "Client" = client
        self._web_host: str = web_host
        self._api_host: str = api_host

    @abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        """将实体序列化为 JSON 可序列化的字典。"""

    def _get(self, path: str, **kwargs) -> Any:
        return self._client.get(path, **kwargs).data

    def _post(self, path: str, **kwargs) -> Any:
        return self._client.post(path, **kwargs).data

    def _put(self, path: str, **kwargs) -> Any:
        return self._client.put(path, **kwargs).data

    def _delete(self, path: str, **kwargs) -> Any:
        return self._client.delete(path, **kwargs).data

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
            items = resp.get("list", []) if isinstance(resp, dict) else resp
            if not items:
                break
            yield from items
            total_pages = resp.get("pages", 1) if isinstance(resp, dict) else 1
            if page >= total_pages:
                break
            page += 1

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        ident = getattr(self, "_path", None) or getattr(self, "_username", None) or "?"
        return f"{cls}('{ident}')"
