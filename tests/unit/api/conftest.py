from unittest.mock import MagicMock

import pytest

from swanlab.api.base import ApiClientContext


@pytest.fixture
def mock_ctx():
    client = MagicMock()
    return ApiClientContext(
        client=client,
        web_host="https://swanlab.cn",
        api_host="https://api.swanlab.cn",
        username="testuser",
        name="Test User",
    )
