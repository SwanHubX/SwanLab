import contextlib
import os
import tempfile
from pathlib import Path
from typing import IO, Any, Optional, Union

from .dir import safe_mkdir


def safe_write(
    target: Union[str, Path, IO[Any]],
    content: Union[str, bytes],
    mode: str = "w",
    atomic: bool = False,
    encoding: Optional[str] = None,
):
    """
    安全地将内容写入目标。
    """
    # 1. 如果传入的是已打开的句柄
    if not isinstance(target, (str, Path)):
        target.write(content)
        target.flush()
        if hasattr(target, "fileno"):
            with contextlib.suppress(OSError):  # 优雅忽略不支持 fsync 的句柄
                os.fsync(target.fileno())
        return

    # 2. 如果传入的是路径
    path = Path(target)

    # 推断编码模式
    is_binary = "b" in mode or isinstance(content, bytes)
    resolved_encoding = None if is_binary else (encoding or "utf-8")
    temp_mode = "wb" if is_binary else "w"

    if atomic:
        # 把算好的确切参数直接喂给底层，底层只管干活
        _atomic_save(path, content, temp_mode, resolved_encoding, is_binary)
    else:
        # 非原子写入也需要保证 NAS 目录就绪！
        safe_mkdir(path.parent)

        with open(path, mode=mode, encoding=resolved_encoding) as f:
            f.write(content)
            f.flush()
            with contextlib.suppress(OSError):
                os.fsync(f.fileno())


def _atomic_save(
    path: Path,
    content: Union[str, bytes],
    temp_mode: str,
    resolved_encoding: Optional[str],
    is_binary: bool,
):
    """
    内部方法：原子级文件保存。
    只负责执行底层 I/O 逻辑，参数已由上层统一解析清洗。
    """
    # 守住 NAS 的防线
    safe_mkdir(path.parent)

    # 在同级目录下创建隐藏的临时文件
    fd, temp_path = tempfile.mkstemp(dir=path.parent, prefix=".tmp_swanlab_", text=not is_binary)

    try:
        # fd 会在 with 块结束时被自动关闭
        with open(fd, temp_mode, encoding=resolved_encoding) as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())

        # 原子替换（POSIX 和现代 Windows 下均安全）
        os.replace(temp_path, path)
    except Exception:
        # 即使删除临时文件失败，也不能掩盖原始异常
        with contextlib.suppress(OSError):
            os.remove(temp_path)
        raise
