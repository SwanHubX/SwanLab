from unittest.mock import MagicMock

from swanlab.core_python.api.experiment import (
    complete_multipart,
    complete_upload,
    prepare_multipart,
    prepare_upload,
)


def test_prepare_upload_wraps_files_payload():
    client = MagicMock()
    files = [{"path": "log/train.log", "size": 12, "md5": "abc123"}]
    client.post.return_value = (
        {"data": {"urls": ["https://upload.example.com/file"]}},
        None,
    )

    result = prepare_upload(client, "exp-123", files)

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/prepare",
        {"files": [{"path": "log/train.log", "size": 12, "md5": "abc123"}]},
    )
    assert result == ["https://upload.example.com/file"]


def test_complete_upload_sends_documented_files_payload():
    client = MagicMock()

    complete_upload(
        client,
        "exp-123",
        [{"path": "log/train.log", "state": "UPLOADED"}],
    )

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/complete",
        {"files": [{"path": "log/train.log", "state": "UPLOADED"}]},
    )


def test_prepare_multipart_sends_documented_fields():
    client = MagicMock()
    client.post.return_value = (
        {
            "data": {
                "files": [
                    {
                        "uploadId": "upload-1",
                        "parts": [
                            {
                                "partNumber": 1,
                                "url": "https://upload.example.com/part-1",
                            }
                        ],
                    }
                ]
            }
        },
        None,
    )

    result = prepare_multipart(
        client,
        "exp-123",
        {
            "path": "checkpoints/model.txt",
            "size": 1024,
            "md5": "deadbeef",
            "count": 3,
            "mimeType": "text/plain",
        },
    )

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/prepare-multipart",
        {
            "files": [
                {
                    "path": "checkpoints/model.txt",
                    "size": 1024,
                    "md5": "deadbeef",
                    "count": 3,
                    "mimeType": "text/plain",
                }
            ]
        },
    )
    assert result == {
        "uploadId": "upload-1",
        "parts": [{"partNumber": 1, "url": "https://upload.example.com/part-1"}],
    }


def test_complete_multipart_sends_upload_id_and_parts():
    client = MagicMock()

    complete_multipart(
        client,
        "exp-123",
        {
            "path": "checkpoints/model.txt",
            "uploadId": "upload-1",
            "parts": [{"partNumber": 1, "etag": "etag-1"}],
        },
    )

    client.post.assert_called_once_with(
        "/experiment/exp-123/files/complete-multipart",
        {
            "files": [
                {
                    "path": "checkpoints/model.txt",
                    "uploadId": "upload-1",
                    "parts": [{"partNumber": 1, "etag": "etag-1"}],
                }
            ]
        },
    )
