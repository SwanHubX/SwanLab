import pytest


@pytest.fixture(autouse=True, scope="function")
def setup_tmp_dir(tmp_path, monkeypatch):
    """
    自动为每个测试用例切换到独立的临时目录
    """
    monkeypatch.chdir(tmp_path)
