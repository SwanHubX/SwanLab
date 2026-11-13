"""四格式真实编解码与逐帧拼贴回归。"""

import subprocess
import weakref
from io import BytesIO
from types import SimpleNamespace

import imageio_ffmpeg
import numpy as np
import pytest
from PIL import Image, ImageSequence
from PIL.GifImagePlugin import GifImageFile

import swanlab.sdk.internal.run.transforms.video as video_module
from swanlab import Video
from swanlab.sdk.internal.run.transforms.video.encoding import MAX_DIMENSION, Frames


def decode(video: Video) -> bytes:
    assert video.is_valid, video.error
    content = video.buffer.getvalue()
    if video.format == "gif":
        with Image.open(BytesIO(content)) as image:
            return b"".join(frame.convert("RGB").tobytes() for frame in ImageSequence.Iterator(image))
    result = subprocess.run(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "error",
            "-i",
            "pipe:0",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgb24",
            "pipe:1",
        ],
        input=content,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
@pytest.mark.parametrize("file_output", [False, True])
@pytest.mark.parametrize("channels", [1, 3, 4])
def test_real_encoding(fmt, file_output, channels, tmp_path, monkeypatch):
    if file_output:
        path = tmp_path / f"encoded.{fmt}"
        path.touch()
        monkeypatch.setattr(video_module, "scratch_path", lambda _: path)
    frames = np.zeros((3, channels, 7, 9), dtype=np.uint8)
    frames[0] = 40
    frames[1] = 120
    frames[2] = 220
    video = Video(frames, format=fmt, fps=4)
    decoded = decode(video)
    width = 32 if fmt == "ogg" else 10 if fmt == "mp4" else 9
    height = 32 if fmt == "ogg" else 8 if fmt == "mp4" else 7
    frame_size = width * height * 3
    assert len(decoded) == frame_size * 3
    pixels = np.frombuffer(decoded[:frame_size], dtype=np.uint8).reshape(height, width, 3)
    assert abs(float(pixels[:7, :9].mean()) - 40) < 15


@pytest.mark.parametrize("batch,rows,cols", [(1, 1, 1), (2, 1, 2), (3, 1, 4), (5, 2, 4), (8, 2, 4)])
@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
def test_batch_layout_and_black_padding(batch, rows, cols, fmt):
    arrays = np.zeros((batch, 2, 3, 4, 6), dtype=np.uint8)
    for index in range(batch):
        arrays[index] = index + 1
    prepared = Frames(arrays, "gif", None)
    frames = list(prepared)
    assert len(frames) == 2
    assert frames[0].shape == (rows * 4, cols * 6, 3)
    for index in range(rows * cols):
        y, x = index // cols * 4, index % cols * 6
        expected = index + 1 if index < batch else 0
        assert np.all(frames[0][y : y + 4, x : x + 6] == expected)
    assert decode(Video(arrays, format=fmt))


@pytest.mark.parametrize(
    "dtype", [bool, complex, object, str, np.int8, np.uint16, np.int32, np.float16, np.float32, np.float64]
)
def test_invalid_array_dtype_is_isolated(dtype):
    video = Video(np.zeros((2, 3, 4, 4), dtype=dtype), format="gif")
    assert not video.is_valid
    assert video.error is not None and video.error.code == "dtype"


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_pixels_are_invalid(value):
    arrays = np.zeros((2, 3, 4, 4))
    arrays[0, 0, 0, 0] = value
    video = Video(arrays, format="gif")
    assert not video.is_valid
    assert video.error is not None and video.error.code == "dtype"


@pytest.mark.parametrize("fps", [0, -1, True, 1.5])
def test_invalid_fps(fps):
    assert not Video(np.zeros((2, 3, 4, 4), dtype=np.uint8), format="gif", fps=fps).is_valid


@pytest.mark.parametrize("protocol", [False, True])
def test_float_pixels_are_rejected_before_encoding(protocol, monkeypatch):
    from swanlab.sdk.internal.run.transforms.video import encoding

    arrays = np.full((1, 3, 4, 4), 0.5, dtype=np.float32)

    class Convertible:
        def numpy(self):
            return arrays

    data = Convertible() if protocol else arrays
    monkeypatch.setattr(encoding, "_encode_gif", lambda *args: pytest.fail("Invalid pixels must not be encoded"))
    video = Video(data, format="gif")
    assert not video.is_valid and video.error is not None
    assert video.error.code == "dtype" and "uint8" in video.error.message
    assert video.transform(step=0, path=None) is None


@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
@pytest.mark.parametrize("shape", [(1, 3, 2, 4097), (1, 3, 4097, 2), (8, 1, 3, 1080, 1920), (8, 1, 3, 2050, 2)])
def test_oversized_canvas_is_rejected_before_encoding(fmt, shape, monkeypatch):
    from swanlab.sdk.internal.run.transforms.video import encoding

    arrays = np.broadcast_to(np.zeros((1,) * len(shape), dtype=np.uint8), shape)
    for name in ("_encode_gif", "_encode_pipe", "_encode_file"):
        monkeypatch.setattr(encoding, name, lambda *args: pytest.fail("Oversized canvas must not be encoded"))
    video = Video(arrays, format=fmt)
    assert not video.is_valid and video.error is not None
    assert video.error.code == "dimensions_limit"
    assert str(MAX_DIMENSION) in video.error.message
    assert video.transform(step=0, path=None) is None


@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
def test_canvas_at_dimension_limit_encodes_without_resizing(fmt):
    arrays = np.zeros((1, 3, 32, MAX_DIMENSION), dtype=np.uint8)
    video = Video(arrays, format=fmt)
    assert len(decode(video)) == 32 * MAX_DIMENSION * 3


@pytest.mark.parametrize("fmt", ["mp4", "ogg"])
def test_padding_and_tiling_respect_dimension_limit(fmt):
    arrays = np.zeros((2, 1, 3, 3, 2047), dtype=np.uint8)
    frames = Frames(arrays, fmt, None)
    assert frames.width <= MAX_DIMENSION and frames.height <= MAX_DIMENSION
    assert next(iter(frames)).shape == (frames.height, frames.width, 3)


def test_numpy_protocol_and_conversion_failure():
    class Convertible:
        def numpy(self):
            return np.zeros((2, 3, 4, 4), dtype=np.uint8)

    class Failing:
        def numpy(self):
            raise RuntimeError("device unavailable")

    assert Video(Convertible(), format="gif").is_valid
    invalid = Video(Failing(), format="gif")
    assert not invalid.is_valid
    assert invalid.error is not None and invalid.error.code == "numpy_conversion"


@pytest.mark.parametrize("fps", [4, 8])
def test_gif_animation_timing(fps):
    frames = np.zeros((3, 3, 4, 4), dtype=np.uint8)
    frames[1] = 100
    frames[2] = 200
    video = Video(frames, format="gif", fps=fps)
    assert video.is_valid
    with Image.open(video.buffer) as image:
        assert isinstance(image, GifImageFile)
        assert image.n_frames == 3
        duration = 0
        for index in range(image.n_frames):
            image.seek(index)
            duration += image.info["duration"]
        assert abs(duration - 3000 / fps) <= 30


def test_dependency_failure_and_interrupt_cleanup(tmp_path, monkeypatch):
    from swanlab.sdk.internal.run.transforms.video import encoding
    from swanlab.sdk.internal.run.transforms.video.source import VideoProcessingError

    path = tmp_path / "owned.gif"
    path.touch()
    monkeypatch.setattr(video_module, "scratch_path", lambda _: path)

    def missing(*args):
        raise VideoProcessingError("dependency_missing", "imageio is missing")

    monkeypatch.setattr(encoding, "_encode_gif", missing)
    invalid = Video(np.zeros((1, 3, 4, 4), dtype=np.uint8), format="gif")
    assert not invalid.is_valid and invalid.error is not None
    assert invalid.error.code == "dependency_missing"
    assert not path.exists()

    def interrupted(*args):
        raise KeyboardInterrupt

    path.touch()
    monkeypatch.setattr(encoding, "_encode_gif", interrupted)
    with pytest.raises(KeyboardInterrupt):
        Video(np.zeros((1, 3, 4, 4), dtype=np.uint8), format="gif")
    assert not path.exists()


@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
def test_existing_file_and_memory_need_no_encoder(fmt, tmp_path, monkeypatch):
    encoded = Video(np.zeros((2, 3, 8, 8), dtype=np.uint8), format=fmt)
    assert encoded.is_valid
    content = encoded.buffer.getvalue()
    path = tmp_path / f"source.{fmt}"
    path.write_bytes(content)

    def no_dependencies(*args):
        pytest.fail("Existing video bytes must not require encoding dependencies")

    monkeypatch.setattr(video_module, "require_modules", no_dependencies)
    file_video = Video(path, fps=0, format="gif")
    memory_video = Video(content, format=fmt)
    assert file_video.is_valid and memory_video.is_valid
    assert file_video.buffer.getvalue() == memory_video.buffer.getvalue() == content


def test_moviepy_error_is_invalid_even_with_decodable_partial_output(tmp_path, monkeypatch):
    from moviepy.video.io import ffmpeg_writer

    arrays = np.zeros((2, 3, 8, 8), dtype=np.uint8)
    encoded = Video(arrays, format="mp4")
    assert encoded.is_valid
    content = encoded.buffer.getvalue()
    path = tmp_path / "owned.mp4"
    path.touch()

    class FailedWriter:
        def __init__(self, filename, *args, **kwargs):
            path.write_bytes(content)
            self.proc = SimpleNamespace(returncode=1)

        def write_frame(self, frame):
            pass

        def close(self):
            pass

    monkeypatch.setattr(ffmpeg_writer, "FFMPEG_VideoWriter", FailedWriter)
    monkeypatch.setattr(video_module, "scratch_path", lambda _: path)
    invalid = Video(arrays, format="mp4")
    assert not invalid.is_valid and invalid.error is not None
    assert invalid.error.code == "encoding"
    assert not path.exists()


@pytest.mark.parametrize("file_output", [False, True])
def test_gif_streaming_retains_at_most_two_rgb_frames(tmp_path, monkeypatch, file_output):
    original = Image.fromarray
    references = []
    live_counts = []

    def fromarray(*args, **kwargs):
        image = original(*args, **kwargs)
        references.append(weakref.ref(image))
        live_counts.append(sum(ref() is not None for ref in references))
        return image

    monkeypatch.setattr(Image, "fromarray", fromarray)
    path = tmp_path / "stream.gif"
    if file_output:
        path.touch()
        monkeypatch.setattr(video_module, "scratch_path", lambda _: path)
    arrays = np.zeros((128, 3, 32, 32), dtype=np.uint8)
    for index in range(len(arrays)):
        arrays[index] = (index % 3) * 80
    video = Video(arrays, format="gif", fps=8)
    assert video.is_valid, video.error
    assert len(references) == 128
    assert max(live_counts) <= 2
    assert len(decode(video)) == 128 * 32 * 32 * 3
    with Image.open(video.buffer) as image:
        assert image.info["loop"] == 0


@pytest.mark.parametrize("batch,rows,cols", [(1, 1, 1), (2, 1, 2), (3, 1, 4), (5, 2, 4), (8, 2, 4)])
@pytest.mark.parametrize("fmt", ["gif", "mp4", "webm", "ogg"])
@pytest.mark.parametrize("file_output", [False, True])
def test_decoded_batch_preserves_positions_black_cells_and_time(
    batch, rows, cols, fmt, file_output, tmp_path, monkeypatch
):
    if file_output:
        path = tmp_path / f"batch.{fmt}"
        path.touch()
        monkeypatch.setattr(video_module, "scratch_path", lambda _: path)
    arrays = np.zeros((batch, 2, 3, 16, 16), dtype=np.uint8)
    for index in range(batch):
        arrays[index, 0] = 40 + index * 20
        arrays[index, 1] = 90 + index * 20
    video = Video(arrays, format=fmt, fps=4)
    height = max(32, rows * 16) if fmt == "ogg" else rows * 16
    width = max(32, cols * 16) if fmt == "ogg" else cols * 16
    decoded = np.frombuffer(decode(video), dtype=np.uint8).reshape(2, height, width, 3)
    for time in range(2):
        for index in range(rows * cols):
            y, x = index // cols * 16, index % cols * 16
            expected = 40 + time * 50 + index * 20 if index < batch else 0
            pixels = decoded[time, y + 3 : y + 13, x + 3 : x + 13]
            assert abs(float(pixels.mean()) - expected) < 15
        if height > rows * 16:
            assert decoded[time, rows * 16 :].mean() < 15
        if width > cols * 16:
            assert decoded[time, :, cols * 16 :].mean() < 15


@pytest.mark.parametrize(
    "source_format,declared,error_code",
    [("mp4", None, "format_required"), ("mp4", "gif", "format_mismatch"), ("gif", "mp4", "format_mismatch")],
)
def test_memory_format_requires_consistent_explicit_container(source_format, declared, error_code):
    encoded = Video(np.zeros((1, 3, 8, 8), dtype=np.uint8), format=source_format)
    assert encoded.is_valid
    invalid = Video(encoded.buffer.getvalue(), format=declared)
    assert not invalid.is_valid and invalid.error is not None
    assert invalid.error.code == error_code


def test_expected_video_warnings_are_emitted(tmp_path, monkeypatch):
    from swanlab.sdk.internal.pkg import console

    messages = []
    monkeypatch.setattr(console, "warning", lambda message: messages.append(message))
    path = tmp_path / "existing.gif"
    path.write_bytes(b"GIF89a-content")
    assert Video(path, fps=8).is_valid
    assert Video(np.zeros((1, 4, 8, 8), dtype=np.uint8), format="gif").is_valid
    assert not Video(np.zeros((1, 3, 8, 8), dtype=np.float32), format="gif").is_valid
    assert any("fps is ignored" in message for message in messages)
    assert any("alpha is not retained" in message for message in messages)
    assert any("dtype uint8" in message for message in messages)
