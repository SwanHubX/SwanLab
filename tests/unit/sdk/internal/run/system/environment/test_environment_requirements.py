"""
@author: cunyue
@file: test_requirements.py
@time: 2026/3/31
@description: Tests for requirements information collection
"""

from subprocess import CalledProcessError, TimeoutExpired
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.run.system.environment import requirements


class TestRequirementsGet:
    """Tests for requirements.get()"""

    def test_get_pixi_success(self):
        """pixi 可用时使用 pixi list"""
        mock_output = "Package         Version  Source\nnumpy           1.26.0   pypi\n"
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=mock_output)
            result = requirements.get()

        assert result == mock_output
        args = mock_run.call_args
        assert args[0][0] == ["pixi", "list"]
        assert args.kwargs.get("timeout") == 5

    def test_get_uv_fallback(self):
        """pixi 不可用时回退到 uv"""
        uv_output = "numpy==1.26.0\npandas==2.1.0\n"
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = [
                MagicMock(returncode=1, stdout=""),
                MagicMock(returncode=0, stdout=uv_output),
            ]
            result = requirements.get()

        assert result == uv_output
        assert mock_run.call_count == 2
        # 验证 uv 参数
        second_call = mock_run.call_args_list[1]
        assert second_call[0][0] == ["uv", "pip", "list", "--format=freeze"]

    def test_get_pip_fallback(self):
        """pixi 和 uv 都不可用时回退到 pip"""
        pip_output = "numpy==1.26.0\npandas==2.1.0\n"
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = [
                MagicMock(returncode=1, stdout=""),  # pixi 失败
                MagicMock(returncode=1, stdout=""),  # uv 失败
                MagicMock(returncode=0, stdout=pip_output),  # pip 成功
            ]
            result = requirements.get()

        assert result == pip_output
        assert mock_run.call_count == 3
        # 验证 pip 参数
        third_call = mock_run.call_args_list[2]
        assert third_call[0][0] == ["pip", "list", "--format=freeze"]
        assert third_call.kwargs.get("timeout") == 15

    @pytest.mark.parametrize(
        "side_effects",
        [
            [
                MagicMock(returncode=1, stdout=""),
                MagicMock(returncode=1, stdout=""),
                CalledProcessError(1, "pip"),
            ],
        ],
    )
    def test_get_returns_none_on_all_failures(self, side_effects):
        """所有工具都失败时返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = side_effects
            result = requirements.get()

        assert result is None

    @pytest.mark.parametrize(
        "exception",
        [
            TimeoutExpired("pip", 15),
            FileNotFoundError("pip not found"),
        ],
    )
    def test_get_returns_none_on_error(self, exception):
        """命令执行异常时返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = exception
            result = requirements.get()

        assert result is None

    @pytest.mark.parametrize(
        "output",
        [
            "",
            "numpy==1.26.0\npandas==2.1.0\nrequests==2.31.0\n",
        ],
    )
    def test_get_output_variations(self, output):
        """测试各种输出情况"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=output)
            result = requirements.get()

        assert result == output
        # 验证通用参数
        kwargs = mock_run.call_args.kwargs
        assert kwargs.get("capture_output") is True
        assert kwargs.get("text") is True
