"""
@author: cunyue
@file: test_conda.py
@time: 2026/3/31
@description: Tests for conda environment information collection
"""

from subprocess import CalledProcessError, TimeoutExpired
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.run.system.environment import conda


class TestCondaGet:
    """Tests for conda.get()"""

    @pytest.mark.parametrize(
        "output,expected",
        [
            ("name: test-env\nchannels:\n  - defaults\n", "name: test-env\nchannels:\n  - defaults\n"),
            ("", ""),
            ("name: myenv\nprefix: /home/user/miniconda3\n", "name: myenv\nprefix: /home/user/miniconda3\n"),
        ],
    )
    def test_get_success(self, output, expected):
        """正常获取 conda 环境信息"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=output)
            result = conda.get()

        assert result == expected
        # 验证调用参数
        kwargs = mock_run.call_args.kwargs
        assert kwargs.get("timeout") == 15
        assert kwargs.get("capture_output") is True
        assert kwargs.get("text") is True

    @pytest.mark.parametrize(
        "exception",
        [
            CalledProcessError(1, "conda"),
            TimeoutExpired("conda", 15),
            FileNotFoundError("conda not found"),
        ],
    )
    def test_get_returns_none_on_error(self, exception):
        """各种错误情况下返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = exception
            result = conda.get()

        assert result is None
