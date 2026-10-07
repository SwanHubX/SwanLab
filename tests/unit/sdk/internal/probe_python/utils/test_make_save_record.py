"""
@author: Nexisato
@file: test_make_save_record.py
@time: 2026/10/7
@description: 测试 probe 内部保存记录的构建（make_save_record）
"""

from pathlib import Path

import pytest

from swanlab.proto.swanlab.save.v1.save_pb2 import SaveType
from swanlab.sdk.internal.probe_python.utils import make_save_record


class TestMakeSaveRecord:
    @pytest.mark.parametrize(
        "name,save_type",
        [
            ("metadata", SaveType.SAVE_TYPE_METADATA),
            ("requirements", SaveType.SAVE_TYPE_REQUIREMENTS),
            ("conda", SaveType.SAVE_TYPE_CONDA),
        ],
    )
    def test_skip_store_without_target(self, name: str, save_type: SaveType):
        record = make_save_record(name, save_type, "内容", None)
        assert record.name == name
        assert record.type == save_type
        assert record.HasField("payload")
        assert record.payload == "内容".encode("utf-8")
        assert record.source_path == ""

    def test_skip_store_payload_inline(self, tmp_path: Path):
        """skip_store：内容进 payload、不落盘、source_path 留空。"""
        target = tmp_path / "metadata.json"
        record = make_save_record("metadata", SaveType.SAVE_TYPE_METADATA, '{"a": 1}', None)
        assert record.name == "metadata"
        assert record.type == SaveType.SAVE_TYPE_METADATA
        assert record.payload == b'{"a": 1}'
        assert record.source_path == ""
        assert not target.exists()

    def test_persist_writes_file(self, tmp_path: Path):
        """默认：写入 target 并以 source_path 引用。"""
        target = tmp_path / "requirements.txt"
        record = make_save_record("requirements", SaveType.SAVE_TYPE_REQUIREMENTS, "pkg==1.0\n", target)
        assert record.name == "requirements"
        assert record.payload == b""
        assert record.source_path == target.absolute().as_posix()
        assert target.read_text(encoding="utf-8") == "pkg==1.0\n"
