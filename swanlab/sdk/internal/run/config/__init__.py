"""
@author: cunyue
@file: __init__.py
@time: 2026/3/14
@description: SwanLab 实验配置模块

生命周期：
  未绑定（bindctx 调用前）：所有写操作仅保留在内存，不触发 IO 和事件
  已绑定（bindctx 调用后）：每次写操作 → 全量覆写 config.yaml → 发出 ConfigEvent
  重置（reset 调用后）    ：清空内存与绑定状态，用于测试隔离或下一次 init

线程安全：
  所有写操作（setitem、delitem、update、pop 等）均在模块级 _lock 内执行。
  读操作（getitem、getattr、iter、len）不加锁（CPython GIL 保证 dict 读原子性）。
  _lock 为普通 Lock（非 RLock），内部方法之间不递归持锁。
"""

import re
import threading
from collections.abc import MutableMapping
from pathlib import Path
from typing import Any, Callable, Optional, Union

from swanlab.proto.swanlab.config.v1.config_pb2 import UpdateType
from swanlab.sdk.internal.bus.events import ConfigEvent

from ._helper import revert_config
from ._parse import parse
from ._writer import write_config

__all__ = [
    "SwanLabConfig",
    "config",
    "create_run_config",
    "create_unbound_run_config",
    "deactivate_run_config",
    "reset",
    "parse",
    "revert_config",
]

# 私有属性名正则：匹配以 _ 开头的属性，拦截用户意外写入
_PRIVATE_RE = re.compile(r"^_")

# 模块级写锁（config 为单例，模块级足够）
_lock = threading.Lock()


class SwanLabConfig(MutableMapping):
    """
    SwanLab 实验配置容器。

    同时支持字典风格（config["lr"] = 0.01）和对象属性风格（config.lr = 0.01）。
    所有写入值均经过 parse() 转换为 JSON 可序列化类型。
    """

    # 类型注解（仅用于 IDE 提示，不创建类属性）
    _config: dict[str, Any]
    _sort: dict[str, int]
    _seq: int
    _file: Optional[Path]
    _emit: Optional[Callable[[ConfigEvent], None]]
    _bound: bool

    def __init__(self) -> None:
        # 直接操作 __dict__ 绕过自定义 __setattr__
        self.__dict__.update(
            {
                "_config": {},  # {key: parsed_value}
                "_sort": {},  # {key: sort_index}
                "_seq": 0,  # 下一个 sort 序号
                "_file": None,  # Path | None
                "_emit": None,  # Callable | None
                "_bound": False,  # 是否已绑定
            }
        )

    # ------------------------------------------------------------------
    # 内部辅助（NOT 线程安全，调用方须持锁）
    # ------------------------------------------------------------------

    def _set_value(self, key: str, value: Any) -> None:
        """解析并写入单个 key，自动维护 sort 序号。"""
        is_new = key not in self._config
        self._config[key] = parse(value)
        if is_new:
            self._sort[key] = self._seq
            self.__dict__["_seq"] = self._seq + 1

    def _flush(self, update_type: UpdateType) -> None:
        """全量写文件并发出 ConfigEvent。"""
        assert self._file is not None and self._emit is not None, "Config not bound"
        write_config(self._file, self._config, self._sort)
        self._emit(ConfigEvent(update=update_type))

    # ------------------------------------------------------------------
    # 绑定 / 重置（线程安全）
    # ------------------------------------------------------------------

    def _bindctx(self, config_file: Path, emit: Callable[[ConfigEvent], None]) -> None:
        """
        绑定运行上下文。将内存中已有的 config 全量 flush 到文件，
        之后的每次写操作均实时同步。可安全重复调用（幂等）。
        """
        with _lock:
            if self._bound:
                return
            self.__dict__.update({"_file": config_file, "_emit": emit, "_bound": True})
            self._flush(UpdateType.UPDATE_TYPE_INIT)

    def _reset(self) -> None:
        """重置为初始状态，用于测试隔离或下一次 init 前的清理。"""
        with _lock:
            self._config.clear()
            self._sort.clear()
            self.__dict__.update({"_seq": 0, "_file": None, "_emit": None, "_bound": False})

    def _snapshot(self) -> tuple[dict, dict, int]:
        """浅拷贝内部状态（config、sort、seq），供 create_run_config 使用。"""
        return dict(self._config), dict(self._sort), self._seq

    def _copy_from(self, source: "SwanLabConfig") -> None:
        """从 source 复制数据，仅在 self 未绑定时调用（无 IO）。"""
        cfg, sort, seq = source._snapshot()
        self.__dict__.update({"_config": cfg, "_sort": sort, "_seq": seq})

    # ------------------------------------------------------------------
    # MutableMapping 接口
    # ------------------------------------------------------------------

    def __setitem__(self, key: str, value: Any) -> None:
        with _lock:
            self._set_value(str(key), value)
            if self._bound:
                self._flush(UpdateType.UPDATE_TYPE_PATCH)

    def __getitem__(self, key: str) -> Any:
        if not isinstance(key, str):
            raise TypeError(f"Key must be a string, got {type(key).__name__}")
        try:
            return self._config[key]
        except KeyError:
            raise KeyError(key)

    def __delitem__(self, key: str) -> None:
        with _lock:
            if key not in self._config:
                raise KeyError(key)
            del self._config[key]
            self._sort.pop(key, None)
            if self._bound:
                self._flush(UpdateType.UPDATE_TYPE_PATCH)

    def __iter__(self):
        return iter(self._config)

    def __len__(self) -> int:
        return len(self._config)

    def __str__(self) -> str:
        return str(self._config)

    # ------------------------------------------------------------------
    # 对象属性风格
    # ------------------------------------------------------------------

    def __setattr__(self, name: str, value: Any) -> None:
        if _PRIVATE_RE.match(name):
            raise AttributeError(f"Attribute '{name}' is private and cannot be set")
        self[name] = value

    def __getattr__(self, name: str) -> Any:
        # 仅在正常属性查找失败时调用；私有字段由 object.__getattribute__ 直接处理
        try:
            return self._config[name]
        except KeyError:
            raise AttributeError(name)

    def __delattr__(self, name: str) -> None:
        if _PRIVATE_RE.match(name):
            raise AttributeError(f"Attribute '{name}' is private and cannot be deleted")
        try:
            del self[name]
        except KeyError:
            raise AttributeError(name)

    # ------------------------------------------------------------------
    # 批量操作（一次 flush）
    # ------------------------------------------------------------------

    def update(self, __m: Optional[Union[MutableMapping, Any]] = None, **kwargs) -> None:  # type: ignore[override]
        """批量更新，所有 key 写完后只触发一次 flush。"""
        with _lock:
            if __m is not None:
                for k, v in parse(__m).items():
                    self._set_value(k, v)
            for k, v in kwargs.items():
                self._set_value(k, v)
            if self._bound:
                self._flush(UpdateType.UPDATE_TYPE_PATCH)

    # ------------------------------------------------------------------
    # 覆盖 MutableMapping 默认实现（避免多余的 flush）
    # ------------------------------------------------------------------

    def get(self, key: str, default: Any = None) -> Any:  # type: ignore[override]
        return self._config.get(key, default)

    def set(self, name: str, value: Any) -> None:
        """显式设置配置项（等价于 config[name] = value）。"""
        self[str(name)] = value

    def pop(self, key: str, *args) -> Any:  # type: ignore[override]
        with _lock:
            if key not in self._config:
                return args[0] if args else None
            value = self._config.pop(key)
            self._sort.pop(key, None)
            if self._bound:
                self._flush(UpdateType.UPDATE_TYPE_PATCH)
            return value

    def clean(self) -> None:
        """清空配置项（不影响绑定状态）。"""
        with _lock:
            self._config.clear()
            self._sort.clear()
            self.__dict__["_seq"] = 0


# ------------------------------------------------------------------
# 全局单例与代理
# ------------------------------------------------------------------


class _ConfigProxy:
    """动态代理：run 活跃时指向 run.config，否则指向 _global_config。"""

    @property
    def _target(self) -> SwanLabConfig:
        return _active_run_config if _active_run_config is not None else _global_config

    def __getitem__(self, key):
        return self._target[key]

    def __setitem__(self, key, value):
        self._target[key] = value

    def __delitem__(self, key):
        del self._target[key]

    def __iter__(self):
        return iter(self._target)

    def __len__(self):
        return len(self._target)

    def __contains__(self, key):
        return key in self._target

    def __str__(self):
        return str(self._target)

    def __getattr__(self, name):
        return getattr(self._target, name)

    def __setattr__(self, name, value):
        setattr(self._target, name, value)

    def __delattr__(self, name):
        delattr(self._target, name)


_global_config = SwanLabConfig()
_active_run_config: Optional[SwanLabConfig] = None

config = _ConfigProxy()


def create_run_config(config_file: Path, emit: Callable) -> SwanLabConfig:
    """从 global config 创建并绑定 per-run config，激活代理。"""
    global _active_run_config
    run_cfg = SwanLabConfig()
    getattr(run_cfg, "_copy_from")(_global_config)
    getattr(run_cfg, "_bindctx")(config_file, emit)
    _active_run_config = run_cfg
    return run_cfg


def create_unbound_run_config() -> SwanLabConfig:
    """disabled 模式：从 global config 创建不绑定文件的 run config，激活代理。"""
    global _active_run_config
    run_cfg = SwanLabConfig()
    getattr(run_cfg, "_copy_from")(_global_config)
    _active_run_config = run_cfg
    return run_cfg


def deactivate_run_config() -> None:
    """run 结束：清理 run config 内存，代理恢复指向 global config。"""
    global _active_run_config
    if _active_run_config is not None:
        getattr(_active_run_config, "_reset")()
        _active_run_config = None


def reset() -> None:
    """重置 global config + 取消激活（仅用于测试隔离）。"""
    global _active_run_config
    getattr(_global_config, "_reset")()
    _active_run_config = None
