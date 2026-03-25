import glob
import os
import pathlib
from typing import List, Union

from .model import SaveFile


def validate_glob_path(glob_path: pathlib.PurePath, base_path: pathlib.PurePath) -> None:
    """
    确保 glob_path 在 base_path 下，避免生成错误的相对路径。
    """
    try:
        glob_path.relative_to(base_path)
    except ValueError as e:
        raise ValueError(f"glob path {glob_path} must be under base path {base_path}") from e


def collect_save_files(
    glob_path: pathlib.PurePath,
    base_path: pathlib.PurePath,
    files_root: Union[str, os.PathLike],
) -> List[SaveFile]:
    """
    根据解析后的 glob/base_path 收集需要保存的文件，并生成 run/files 下的目标路径。
    """
    relative_glob = glob_path.relative_to(base_path)
    source_glob = str(base_path / relative_glob)
    matched = sorted(glob.glob(source_glob, recursive=True))
    files_root = pathlib.Path(files_root)

    files: List[SaveFile] = []
    seen_paths = set()
    seen_names = {}

    for path in matched:
        src = pathlib.Path(path).absolute()
        if str(src) in seen_paths:
            continue
        if not src.is_file():
            continue

        rel = pathlib.Path(*src.parts[len(base_path.parts) :])
        name = rel.as_posix()
        existed = seen_names.get(name)
        if existed is not None and existed != str(src):
            raise ValueError(f'Multiple files resolve to the same save key "{name}".')

        seen_paths.add(str(src))
        seen_names[name] = str(src)
        files.append(
            SaveFile(
                source_path=str(src),
                name=name,
                target_path=str((files_root / rel).absolute()),
            )
        )

    return files


__all__ = [
    "SaveFile",
    "validate_glob_path",
    "collect_save_files",
]
