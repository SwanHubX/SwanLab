"""
@author: likang
@file: test_login.py
@time: 2026/11/4
@description: swanlab login CLI 命令的错误退出行为
"""

import requests
import responses
from click.testing import CliRunner

from swanlab.cli.auth import login as login_command


@responses.activate
def test_network_error_exits_without_traceback():
    """CLI 命令：连接失败应给出可读提示并以非 0 退出，而不是抛 traceback。"""
    responses.add(
        responses.POST,
        "https://api.error.test/api/login/api_key",
        body=requests.ConnectionError("boom"),
    )

    result = CliRunner().invoke(login_command, ["--api-key", "test-key", "--host", "https://api.error.test"])

    assert result.exit_code == 1
    assert "Traceback" not in result.output
    assert "Traceback" not in getattr(result, "stderr", "")
