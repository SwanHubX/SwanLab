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
def test_network_error_exits_with_non_zero_code():
    """CLI 命令：连接失败应给出可读提示并以非 0 退出，而不是让异常逃逸。"""
    responses.add(
        responses.POST,
        "https://api.error.test/api/login/api_key",
        body=requests.ConnectionError("boom"),
    )

    result = CliRunner().invoke(login_command, ["--api-key", "test-key", "--host", "https://api.error.test"])

    assert result.exit_code == 1
    # 命令层捕获后主动退出，而不是未处理异常逃逸
    assert isinstance(result.exception, SystemExit)
