from unittest.mock import MagicMock

from swanlab.core_python.api.experiment import (
    complete_multipart,
    complete_upload,
    prepare_multipart,
    prepare_upload,
)
from swanlab.core_python.save import SaveFilePayload, SaveFileState


def test_prepare_upload_wraps_files_payload():
    client = MagicMock()
    files = [SaveFilePayload(name="log/train.log", size=12, md5="abc123")]
    client.post.return_value = ({"data": {"files": [{"uploadUrl": "https://upload.example.com/file"}]}}, None)

    result = prepare_upload(client, "exp-123", files)

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/prepare",
        {"files": [{"name": "log/train.log", "size": 12, "md5": "abc123"}]},
    )
    assert result == [{"uploadUrl": "https://upload.example.com/file"}]


def test_complete_upload_sends_documented_files_payload():
    client = MagicMock()

    complete_upload(client, "exp-123", [SaveFilePayload(name="log/train.log", state=SaveFileState.UPLOADED)])

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/complete",
        {"files": [{"name": "log/train.log", "state": "UPLOADED"}]},
    )


def test_prepare_multipart_sends_documented_fields():
    client = MagicMock()
    client.post.return_value = (
        {"data": {"uploadId": "upload-1", "uploadUrls": ["https://upload.example.com/part-1"]}},
        None,
    )

    result = prepare_multipart(
        client,
        "exp-123",
        SaveFilePayload(
            name="checkpoints/model.txt",
            size=1024,
            md5="deadbeef",
            count=3,
            mime_type="text/plain",
        ),
    )

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/prepare-multipart",
        {
            "files": [
                {
                    "name": "checkpoints/model.txt",
                    "size": 1024,
                    "md5": "deadbeef",
                    "count": 3,
                    "mimeType": "text/plain",
                }
            ]
        },
    )
    assert result == {"uploadId": "upload-1", "uploadUrls": ["https://upload.example.com/part-1"]}


def test_complete_multipart_sends_upload_id_and_state():
    client = MagicMock()

    complete_multipart(
        client,
        "exp-123",
        SaveFilePayload(
            name="checkpoints/model.txt",
            upload_id="upload-1",
            state=SaveFileState.UPLOADED,
        ),
    )

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/complete-multipart",
        {"files": [{"name": "checkpoints/model.txt", "uploadId": "upload-1", "state": "UPLOADED"}]},
    )
