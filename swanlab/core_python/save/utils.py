import glob
import hashlib
import mimetypes
import os
import pathlib
from typing import Iterable, List, Optional, Union

from .model import FileSignature, SaveFileModel

MIME_TYPE_DEFAULT: str = "application/octet-stream"


def _expand_matched_path(path: Union[str, os.PathLike]) -> Iterable[pathlib.Path]:
    src = pathlib.Path(path).absolute()
    if src.is_file():
        return [src]
    if src.is_dir():
        return sorted(child.absolute() for child in src.rglob("*") if child.is_file())
    return []


def validate_glob_path(glob_path: pathlib.PurePath, base_path: pathlib.PurePath) -> None:
    """确保 glob_path 在 base_path 下"""
    try:
        glob_path.relative_to(base_path)
    except ValueError as e:
        raise ValueError(f"glob path {glob_path} must be under base path {base_path}") from e


def collect_save_files(
    glob_path: pathlib.PurePath,
    base_path: pathlib.PurePath,
    files_root: Union[str, os.PathLike],
) -> List[SaveFileModel]:
    """收集需要保存的文件"""
    relative_glob = glob_path.relative_to(base_path)
    matched = sorted(glob.glob(str(base_path / relative_glob), recursive=True))
    files_root = pathlib.Path(files_root)

    files: List[SaveFileModel] = []
    seen_paths = set()
    seen_names = {}

    for path in matched:
        for src in _expand_matched_path(path):
            if str(src) in seen_paths:
                continue

            rel = pathlib.Path(*src.parts[len(base_path.parts) :])
            name = rel.as_posix()
            if name in seen_names and seen_names[name] != str(src):
                raise ValueError(f'Multiple files resolve to the same save key "{name}".')

            seen_paths.add(str(src))
            seen_names[name] = str(src)
            files.append(
                SaveFileModel(
                    source_path=str(src),
                    name=name,
                    target_path=str((files_root / rel).absolute()),
                )
            )

    return files


def compute_md5(path: str, chunk_size: int = 1024 * 1024) -> str:
    """计算文件 MD5"""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def guess_mime_type(path: str) -> str:
    """推测文件 MIME 类型"""
    mime_type, _ = mimetypes.guess_type(path)
    return mime_type or MIME_TYPE_DEFAULT


def file_signature(path: str) -> Optional[FileSignature]:
    """获取文件签名 (mtime_ns, size)"""
    try:
        stat = os.stat(path)
        return (stat.st_mtime_ns, stat.st_size)
    except OSError:
        return None


__all__ = [
    "validate_glob_path",
    "collect_save_files",
    "compute_md5",
    "guess_mime_type",
    "file_signature",
]
