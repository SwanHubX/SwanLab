"""
@author: cunyue
@file: test_media_transform.py
@time: 2026/3/15
@description: 所有 TransformMedia 子类的 transform 通用契约测试。

新增媒体类型时，只需在 MEDIA_FACTORIES 中追加一项；
若忘记注册，test_all_media_types_have_factory 会报错提示。
"""

import inspect
import re
from pathlib import Path

import numpy as np
import pytest

import swanlab.sdk.internal.run.transforms  # noqa: F401  # type: ignore — 触发所有子类注册
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.run.transforms.audio import Audio
from swanlab.sdk.internal.run.transforms.image import Image
from swanlab.sdk.internal.run.transforms.text import Text

# 注册表：TransformMedia 子类 → 无参工厂（每次调用返回内容相同的新实例）
MEDIA_FACTORIES = {
    Audio: lambda: Audio(np.zeros((1, 4410), dtype=np.float32), sample_rate=44100),
    Image: lambda: Image(np.zeros((10, 10, 3), dtype=np.uint8)),
    Text: lambda: Text(content="hello world"),
}


def _discover_concrete_subclasses(cls):
    """递归发现所有非抽象子类"""
    result = []
    for sub in cls.__subclasses__():
        if not inspect.isabstract(sub):
            result.append(sub)
        result.extend(_discover_concrete_subclasses(sub))
    return result


def test_all_media_types_have_factory():
    """防护测试：新增 TransformMedia 子类但忘记在 MEDIA_FACTORIES 注册工厂时，此测试报错"""
    all_concrete = _discover_concrete_subclasses(TransformMedia)
    unregistered = [cls for cls in all_concrete if cls not in MEDIA_FACTORIES]
    assert not unregistered, "以下 TransformMedia 子类缺少测试工厂，请在 MEDIA_FACTORIES 中注册:\n" + "\n".join(
        f"  - {cls.__qualname__}" for cls in unregistered
    )


@pytest.mark.parametrize(
    "factory",
    [pytest.param(f, id=cls.__name__) for cls, f in MEDIA_FACTORIES.items()],
)
class TestMediaTransformContract:
    """
    transform 通用契约：
      1. 文件写入磁盘
      2. 文件名格式为 {step:03d}-{8位小写十六进制}.{ext}
      3. step 不同 → 文件名不同，且前缀正确反映 step
      4. 相同内容 + 相同 step → 相同文件名（内容寻址）
    """

    def test_file_written_to_disk(self, factory, tmp_path: Path):
        item = factory().transform(step=1, path=tmp_path)
        assert (tmp_path / item.filename).exists()

    def test_filename_format(self, factory, tmp_path: Path):
        item = factory().transform(step=7, path=tmp_path)
        assert re.match(r"^007-[0-9a-f]{8}\.", item.filename), f"Unexpected filename: {item.filename!r}"

    def test_step_in_filename(self, factory, tmp_path: Path):
        m = factory()
        item1 = m.transform(step=1, path=tmp_path)
        item2 = m.transform(step=2, path=tmp_path)
        assert item1.filename.startswith("001-")
        assert item2.filename.startswith("002-")
        assert item1.filename != item2.filename

    def test_same_content_same_filename(self, factory, tmp_path: Path):
        item1 = factory().transform(step=1, path=tmp_path)
        item2 = factory().transform(step=1, path=tmp_path)
        assert item1.filename == item2.filename
