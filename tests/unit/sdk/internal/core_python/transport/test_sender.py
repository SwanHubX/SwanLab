from pathlib import Path
from unittest.mock import MagicMock, patch

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.settings.core import SaveSettings


def _make_sender(tmp_path: Path) -> HttpRecordSender:
    config = MagicMock()
    config.settings.core.save = SaveSettings()
    ctx = MagicMock(spec=RunContext)
    ctx.config = config
    return HttpRecordSender(
        run_dir=tmp_path,
        username="alice",
        project="demo",
        project_id="project-id",
        experiment_id="experiment-id",
        ctx=ctx,
    )


def _make_save_record(source: Path, name: str = "checkpoints/model.txt") -> Record:
    return Record(save=SaveRecord(name=name, source_path=str(source), target_path=str(source)))


def test_upload_save_delegates_small_file_s3_put_to_upload_resources(tmp_path: Path):
    source = tmp_path / "model.txt"
    source.write_text("weights", encoding="utf-8")
    sender = _make_sender(tmp_path)
    session = MagicMock()

    def _fake_upload(session, *, resources):
        return {r["source_path"] for r in resources}

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/model"]},
        ) as mock_prepare,
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=session),
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.upload_saves",
            side_effect=_fake_upload,
        ) as mock_upload_saves,
    ):
        sender.upload_save([_make_save_record(source)])

    mock_prepare.assert_called_once()
    mock_upload_saves.assert_called_once_with(
        session,
        resources=[
            {"url": "https://s3.test/model", "source_path": str(source), "content_type": "text/plain"},
        ],
    )
    mock_complete.assert_called_once_with(
        "experiment-id", files=[{"path": "checkpoints/model.txt", "state": "UPLOADED"}]
    )


def test_upload_save_skips_complete_when_s3_upload_fails(tmp_path: Path):
    source = tmp_path / "failed.bin"
    source.write_bytes(b"failed")
    sender = _make_sender(tmp_path)

    with (
        patch(
            "swanlab.sdk.internal.core_python.transport.sender.prepare_save_files",
            return_value={"urls": ["https://s3.test/failed"]},
        ),
        patch("swanlab.sdk.internal.core_python.transport.sender.complete_save_files") as mock_complete,
        patch("swanlab.sdk.internal.core_python.transport.sender.client.session.create", return_value=MagicMock()),
        patch("swanlab.sdk.internal.core_python.transport.sender.upload_saves", side_effect=RuntimeError("s3 failed")),
        patch("swanlab.sdk.internal.pkg.safe.console.trace"),
    ):
        sender.upload_save([_make_save_record(source, name="failed.bin")])

    # 上传失败时仍应调用 complete，但标记为 FAILED
    mock_complete.assert_called_once_with("experiment-id", files=[{"path": "failed.bin", "state": "FAILED"}])
