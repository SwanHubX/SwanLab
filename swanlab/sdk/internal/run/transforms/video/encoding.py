"""逐帧准备视频，使用 MoviePy 文件 writer 或 FFmpeg 内存管道编码。"""

from __future__ import annotations

import importlib.util
import subprocess
import threading
from contextlib import nullcontext
from io import BytesIO
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Optional, cast

from swanlab.sdk.internal.pkg import console

from .source import CHUNK_SIZE, MAX_SIZE, VideoProcessingError, VideoSource

if TYPE_CHECKING:
    import numpy as np

CODECS = {"mp4": "libx264", "webm": "libvpx", "ogg": "libtheora"}
MAX_DIMENSION = 4096


def require_modules(*names: str) -> None:
    for name in names:
        if importlib.util.find_spec(name) is None:
            raise VideoProcessingError("dependency_missing", f"Video encoding requires {name}; install swanlab[media]")


class Frames:
    def __init__(self, data: np.ndarray, fmt: str, fps: Optional[int]):
        import numpy as np

        if not isinstance(data, np.ndarray):
            raise VideoProcessingError("numpy_conversion", "numpy() must return a NumPy ndarray")
        if data.ndim not in (4, 5) or any(n <= 0 for n in data.shape):
            raise VideoProcessingError("shape", "Video requires positive TCHW or BTCHW dimensions")
        if data.dtype != np.uint8:
            raise VideoProcessingError("dtype", "Video pixels must have dtype uint8; convert the input explicitly")
        if fps is not None and (isinstance(fps, bool) or not isinstance(fps, Integral) or fps <= 0):
            raise VideoProcessingError("argument", "Video fps must be a positive non-boolean integer")
        self.fps = int(fps) if fps is not None else 4
        self.data = data.reshape((1, *data.shape)) if data.ndim == 4 else data
        self.batch, self.count, channels, height, width = tuple(int(n) for n in self.data.shape)
        if channels not in (1, 3, 4):
            raise VideoProcessingError("shape", "Video channels must be 1, 3 or 4")
        if channels == 4:
            console.warning("Video RGBA frames are encoded as RGB; alpha is not retained")
        self.rows = 2 ** ((self.batch.bit_length() - 1) // 2)
        padded_batch = 1 << (self.batch - 1).bit_length()
        self.cols = padded_batch // self.rows
        self.content_width, self.content_height = width * self.cols, height * self.rows
        self.width, self.height = self.content_width, self.content_height
        if fmt == "mp4":
            self.width += self.width % 2
            self.height += self.height % 2
        elif fmt == "ogg":
            self.width = max(32, (self.width + 15) // 16 * 16)
            self.height = max(32, (self.height + 15) // 16 * 16)
        if self.width > MAX_DIMENSION or self.height > MAX_DIMENSION:
            raise VideoProcessingError(
                "dimensions_limit",
                f"Video canvas {self.width}x{self.height} exceeds the {MAX_DIMENSION}-pixel limit per side",
            )

    def __iter__(self) -> Iterator[np.ndarray]:
        import numpy as np

        _, _, channels, height, width = tuple(int(n) for n in self.data.shape)
        for index in range(self.count):
            canvas = np.zeros((self.height, self.width, 3), dtype=np.uint8)
            for batch in range(self.batch):
                frame = self.data[batch, index]
                rgb = np.transpose(frame, (1, 2, 0))
                if channels == 1:
                    rgb = np.repeat(rgb, 3, axis=2)
                y, x = batch // self.cols * height, batch % self.cols * width
                canvas[y : y + height, x : x + width] = rgb[:, :, :3]
            yield canvas


def _binary() -> str:
    require_modules("imageio_ffmpeg")
    import imageio_ffmpeg

    # 不调用下载接口；get_ffmpeg_exe 仅查找已有二进制。
    return imageio_ffmpeg.get_ffmpeg_exe()


def _popen_options() -> dict:
    return {"creationflags": subprocess.CREATE_NO_WINDOW} if hasattr(subprocess, "CREATE_NO_WINDOW") else {}


def _encode_pipe(frames: Frames, fmt: str) -> BytesIO:
    command = [
        _binary(),
        "-hide_banner",
        "-loglevel",
        "error",
        "-f",
        "rawvideo",
        "-pix_fmt",
        "rgb24",
        "-s",
        f"{frames.width}x{frames.height}",
        "-r",
        str(frames.fps),
        "-i",
        "pipe:0",
        "-an",
        "-c:v",
        CODECS[fmt],
        "-pix_fmt",
        "yuv420p",
    ]
    if fmt == "mp4":
        command += ["-movflags", "frag_keyframe+empty_moov"]
    command += ["-f", {"mp4": "mp4", "webm": "webm", "ogg": "ogg"}[fmt], "pipe:1"]
    output, diagnostics = BytesIO(), bytearray()
    exceeded = threading.Event()
    errors: list[Exception] = []
    process = cast(
        subprocess.Popen[bytes],
        subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **_popen_options()
        ),
    )
    assert process.stdin is not None and process.stdout is not None and process.stderr is not None

    def drain_output() -> None:
        assert process.stdout is not None
        try:
            while True:
                chunk = process.stdout.read(CHUNK_SIZE)
                if not chunk:
                    break
                if output.tell() + len(chunk) > MAX_SIZE:
                    exceeded.set()
                    process.terminate()
                    break
                output.write(chunk)
        except Exception as error:
            errors.append(error)
            process.terminate()

    def drain_errors() -> None:
        assert process.stderr is not None
        try:
            while True:
                chunk = process.stderr.read(4096)
                if not chunk:
                    break
                diagnostics.extend(chunk)
                del diagnostics[:-8192]
        except Exception as error:
            errors.append(error)

    readers = [threading.Thread(target=drain_output), threading.Thread(target=drain_errors)]
    for thread in readers:
        thread.start()
    try:
        for frame in frames:
            process.stdin.write(frame.tobytes())
        process.stdin.close()
        process.wait()
    except BrokenPipeError:
        pass
    finally:
        if process.poll() is None:
            process.terminate()
        try:
            process.stdin.close()
        except OSError:
            pass
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        for thread in readers:
            thread.join()
        process.stdout.close()
        process.stderr.close()
    if exceeded.is_set():
        raise VideoProcessingError("size_limit", f"Encoded video exceeds {MAX_SIZE} bytes")
    if process.returncode or errors:
        reason = diagnostics.decode("utf-8", errors="replace").strip()[-1000:]
        raise VideoProcessingError("encoding", f"FFmpeg {fmt} encoder failed: {reason}")
    output.seek(0)
    return output


def _encode_file(frames: Frames, fmt: str, path: Path) -> None:
    require_modules("moviepy", "imageio_ffmpeg")
    from moviepy.video.io.ffmpeg_writer import FFMPEG_VideoWriter

    # 直接使用 MoviePy writer，逐帧传入；无需持有完整 ImageSequenceClip。
    writer = FFMPEG_VideoWriter(
        str(path),
        (frames.width, frames.height),
        frames.fps,
        codec=CODECS[fmt],
        ffmpeg_params=["-pix_fmt", "yuv420p"],
    )
    process = writer.proc
    if process is None:
        writer.close()
        raise VideoProcessingError("encoding", "MoviePy did not start an encoder process")
    try:
        for frame in frames:
            writer.write_frame(frame)
            if path.stat().st_size > MAX_SIZE:
                raise VideoProcessingError("size_limit", f"Encoded video exceeds {MAX_SIZE} bytes")
    finally:
        writer.close()
    if process.returncode:
        raise VideoProcessingError("encoding", f"MoviePy {fmt} encoder exited with an error")


def _encode_gif(frames: Frames, path: Optional[Path]) -> Optional[BytesIO]:
    require_modules("PIL")
    from PIL import GifImagePlugin, Image

    memory = BytesIO() if path is None else None
    with path.open("wb") if path is not None else nullcontext(memory) as output:
        assert output is not None
        size = 0

        def write(block: bytes) -> None:
            nonlocal size
            size += len(block)
            if size > MAX_SIZE:
                raise VideoProcessingError("size_limit", f"Encoded video exceeds {MAX_SIZE} bytes")
            output.write(block)

        for index, frame in enumerate(frames):
            # 每帧单独量化并写入，及时关闭 PIL 图像；不保留整段未压缩帧。
            with Image.fromarray(frame) as rgb, rgb.quantize(colors=256, method=0) as image:
                if index == 0:
                    header, _ = GifImagePlugin.getheader(image, info={"loop": 0, "optimize": False})
                    for block in header:
                        write(block)
                blocks = GifImagePlugin.getdata(
                    image,
                    duration=1000 / frames.fps,
                    disposal=2,
                    include_color_table=index > 0,
                )
                for block in blocks:
                    write(block)
        write(b";")
    if memory is not None:
        memory.seek(0)
    return memory


def _verify(source: VideoSource, frames: Frames) -> None:
    if source.format == "gif":
        from PIL import Image

        with source.reader() as stream, Image.open(stream) as image:
            image.load()
            if image.size != (frames.width, frames.height):
                raise VideoProcessingError("encoded_file_invalid", "Encoded GIF dimensions do not match its canvas")
        return
    command = [_binary(), "-hide_banner", "-loglevel", "error", "-i"]
    memory = source.memory
    command += [
        str(source.path) if memory is None else "pipe:0",
        "-frames:v",
        "1",
        "-f",
        "rawvideo",
        "-pix_fmt",
        "rgb24",
        "pipe:1",
    ]
    view = memory.getbuffer() if memory is not None else None
    try:
        result = subprocess.run(command, input=view, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **_popen_options())
    finally:
        if view is not None:
            view.release()
    if result.returncode or len(result.stdout) != frames.width * frames.height * 3:
        raise VideoProcessingError("encoded_file_invalid", f"Encoded {source.format} video cannot be decoded")


def encode(data: np.ndarray, fmt: str, fps: Optional[int], path: Optional[Path]) -> VideoSource:
    source = VideoSource(fmt, path=path)
    try:
        frames = Frames(data, fmt, fps)
        if fmt == "gif":
            source.memory = _encode_gif(frames, path)
        elif path is None:
            source.memory = _encode_pipe(frames, fmt)
        else:
            _encode_file(frames, fmt, path)
        source.refresh()
        _verify(source, frames)
        return source
    except BaseException as error:
        if source._cleanup is not None:
            source._cleanup()
        if isinstance(error, Exception) and not isinstance(error, VideoProcessingError):
            raise VideoProcessingError("encoding", f"Video encoding failed: {type(error).__name__}: {error}") from error
        raise
