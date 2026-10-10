import pytest
import responses

from swanlab.sdk.internal.core_python import client


@pytest.fixture
def core_auth():
    """Mock verification only; Core still owns the real runtime Client."""
    host = "https://core.example.invalid"
    with responses.RequestsMock(assert_all_requests_are_fired=False) as requests:
        requests.add(
            responses.GET,
            f"{host}/api/auth/verify",
            json={"uid": 1, "username": "test-user", "createdAt": "2026-01-01T00:00:00Z"},
        )
        yield host
        client.reset()
