"""core_go 测试共享夹具。"""

import shutil
import sys
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def short_tmp_path(tmp_path):
    """返回足够短的临时目录，避免 Unix socket 路径超限。

    POSIX 上 Unix domain socket 路径通常限制在 104 字节以内，pytest 的
    ``tmp_path`` 可能超过该限制；Windows 仍使用 ``tmp_path``。
    """
    if sys.platform == "win32":
        yield tmp_path
        return

    root = Path(tempfile.mkdtemp(prefix="swl-t-", dir="/tmp"))
    try:
        yield root
    finally:
        shutil.rmtree(root, ignore_errors=True)
