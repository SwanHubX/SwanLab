"""视频数据源：受管快照、共享内存、分段校验与哈希。"""

from __future__ import annotations

import hashlib
import weakref
from contextlib import contextmanager, nullcontext
from io import BytesIO
from pathlib import Path
from typing import BinaryIO, Iterator, Optional
from uuid import uuid4

from swanlab.sdk.internal.pkg import console

CHUNK_SIZE = 1024 * 1024
MAX_SIZE = 100 * 1024 * 1024
SUPPORTED_FORMATS = ("gif", "mp4", "webm", "ogg")


class VideoProcessingError(ValueError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def check_size(size: int) -> None:
    if size > MAX_SIZE:
        raise VideoProcessingError("size_limit", f"Video exceeds {MAX_SIZE} bytes")
    if size <= 0:
        raise VideoProcessingError("format_mismatch", "Video content is empty")


def _vint(data: bytes, offset: int, *, identifier: bool = False) -> tuple[int, int]:
    if offset >= len(data) or data[offset] == 0:
        raise VideoProcessingError("format_mismatch", "Invalid WebM EBML header")
    length = 1
    mask = 0x80
    while not data[offset] & mask:
        mask >>= 1
        length += 1
    if length > 8 or offset + length > len(data):
        raise VideoProcessingError("format_mismatch", "Truncated WebM EBML header")
    value = int.from_bytes(data[offset : offset + length], "big")
    if not identifier:
        value &= (1 << (7 * length)) - 1
    return value, offset + length


def validate_header(header: bytes, fmt: str) -> None:
    valid = False
    if fmt == "gif":
        valid = header.startswith((b"GIF87a", b"GIF89a"))
    elif fmt == "ogg":
        valid = header.startswith(b"OggS")
    elif fmt == "mp4":
        offset = 0
        while offset + 8 <= len(header):
            size = int.from_bytes(header[offset : offset + 4], "big")
            kind = header[offset + 4 : offset + 8]
            prefix = 8
            if size == 1:
                if offset + 16 > len(header):
                    break
                size = int.from_bytes(header[offset + 8 : offset + 16], "big")
                prefix = 16
            if size < prefix or offset + size > len(header):
                break
            if kind == b"ftyp" and size >= prefix + 8:
                brands = [header[offset + prefix : offset + prefix + 4]]
                brands.extend(header[i : i + 4] for i in range(offset + prefix + 8, offset + size - 3, 4))
                valid = any(
                    b.startswith(b"iso") or b in (b"mp41", b"mp42", b"avc1", b"M4V ", b"dash", b"MSNV") for b in brands
                )
                break
            offset += size
    elif fmt == "webm" and header.startswith(b"\x1a\x45\xdf\xa3"):
        length, offset = _vint(header, 4)
        end = min(offset + length, len(header))
        while offset < end:
            kind, offset = _vint(header, offset, identifier=True)
            size, offset = _vint(header, offset)
            if offset + size > end:
                break
            if kind == 0x4282:
                valid = header[offset : offset + size] == b"webm"
                break
            offset += size
    if not valid:
        raise VideoProcessingError("format_mismatch", f"Video header does not match {fmt}")


def _remove_scratch(path: Path) -> None:
    # 只接收当前对象创建的随机 scratch 路径，不枚举目录。
    try:
        path.unlink(missing_ok=True)
    except OSError:
        console.warning("Unable to remove an SDK video scratch file")


def scratch_path(fmt: str) -> Optional[Path]:
    from swanlab.sdk.internal.run import get_run, has_run

    if not has_run():
        return None
    ctx = get_run()._ctx
    if ctx.config.settings.mode == "disabled" or ctx.config.settings.core.skip_store:
        return None
    directory = ctx.media_dir / "video" / ".scratch"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{uuid4().hex}.{fmt}"
    # 排他创建，后续 writer 仅能覆盖本对象拥有的文件。
    with path.open("xb"):
        pass
    return path


class VideoSource:
    def __init__(self, fmt: str, *, path: Optional[Path] = None, buffer: Optional[BytesIO] = None):
        self.format = fmt
        self.path = path
        self.memory = buffer
        self.size = 0
        self.sha256 = ""
        self._cleanup = weakref.finalize(self, _remove_scratch, path) if path is not None else None

    @contextmanager
    def reader(self) -> Iterator[BinaryIO]:
        if self.memory is not None:
            position = self.memory.tell()
            self.memory.seek(0)
            try:
                yield self.memory
            finally:
                self.memory.seek(position)
        elif self.path is not None:
            with self.path.open("rb") as stream:
                yield stream
        else:
            raise VideoProcessingError("file_read", "Video data source is unavailable")

    def refresh(self) -> None:
        digest = hashlib.sha256()
        size = 0
        with self.reader() as stream:
            header = stream.read(65536)
            validate_header(header, self.format)
            digest.update(header)
            size += len(header)
            while True:
                chunk = stream.read(CHUNK_SIZE)
                if not chunk:
                    break
                size += len(chunk)
                check_size(size)
                digest.update(chunk)
        check_size(size)
        sha256 = digest.hexdigest()
        if self.sha256 and sha256 != self.sha256:
            console.warning("Video content changed; its size and hash were refreshed")
        self.size, self.sha256 = size, sha256

    @property
    def buffer(self) -> BytesIO:
        if self.memory is None:
            result = BytesIO()
            with self.reader() as stream:
                total = 0
                while True:
                    chunk = stream.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    total += len(chunk)
                    check_size(total)
                    result.write(chunk)
            result.seek(0)
            self.memory = result
        return self.memory

    def store(self, target: Path) -> None:
        if self.path == target and self.memory is None:
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        # 已有同 hash 文件也需验证，避免覆盖不属于本对象的不同内容。
        if target.exists():
            digest = hashlib.sha256()
            with target.open("rb") as existing:
                for chunk in iter(lambda: existing.read(CHUNK_SIZE), b""):
                    digest.update(chunk)
            if digest.hexdigest() != self.sha256:
                raise VideoProcessingError("file_read", "Video destination has unexpected content")
        else:
            with self.reader() as source, target.open("xb") as output:
                total = 0
                digest = hashlib.sha256()
                for chunk in iter(lambda: source.read(CHUNK_SIZE), b""):
                    total += len(chunk)
                    check_size(total)
                    digest.update(chunk)
                    output.write(chunk)
                if total != self.size or digest.hexdigest() != self.sha256:
                    raise VideoProcessingError("content_changed", "Video source changed while being saved")
        if self._cleanup is not None and self._cleanup.alive:
            self._cleanup()
        self.path = target


def snapshot(stream: BinaryIO, fmt: str, path: Optional[Path]) -> VideoSource:
    source = VideoSource(fmt, path=path, buffer=BytesIO() if path is None else None)
    try:
        digest = hashlib.sha256()
        header = bytearray()
        with path.open("wb") if path is not None else nullcontext(source.memory) as output:
            assert output is not None
            size = 0
            while True:
                chunk = stream.read(min(CHUNK_SIZE, MAX_SIZE + 1 - size))
                if not chunk:
                    break
                size += len(chunk)
                check_size(size)
                if len(header) < 65536:
                    header.extend(chunk[: 65536 - len(header)])
                digest.update(chunk)
                output.write(chunk)
        check_size(size)
        validate_header(bytes(header), fmt)
        source.size = size
        source.sha256 = digest.hexdigest()
        if source.memory is not None:
            source.memory.seek(0)
        return source
    except BaseException:
        if source._cleanup is not None:
            source._cleanup()
        raise
