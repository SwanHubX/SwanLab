"""
@file: conftest.py
@description: Pytest configuration for converter tests
"""

import os
import pytest


@pytest.fixture(scope="session", autouse=True)
def setup_pytest_env():
    """Set PYTEST_VERSION to allow reset_run_store() calls in teardown"""
    os.environ['PYTEST_VERSION'] = '1'
    yield
    if 'PYTEST_VERSION' in os.environ:
        del os.environ['PYTEST_VERSION']
