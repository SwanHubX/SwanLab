"""
@author: cunyue
@file: test_service.py
@time: 2025/12/11 19:42
@description: 测试服务相关函数
"""

from io import BytesIO
from unittest.mock import MagicMock, patch, call

import pytest
import requests

from swanlab.core_python.api.service import upload_file


@pytest.fixture
def mock_buffer():
    return BytesIO(b"test data content")


@pytest.fixture
def target_url():
    return "https://cos.example.com/upload"


class TestUpload:

    @patch('requests.Session')
    def test_upload_success_first_try(self, mock_session_cls, mock_buffer, target_url):
        """测试场景1：一次直接上传成功"""

        # 1. 配置 Mock
        mock_session_instance = mock_session_cls.return_value
        mock_adapter = mock_session_instance.__enter__.return_value

        # 模拟 put 返回 200 OK
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_adapter.put.return_value = mock_response

        # 2. 执行函数
        upload_file(url=target_url, buffer=mock_buffer)

        # 3. 验证
        # 验证 put 是否被调用了一次
        mock_adapter.put.assert_called_once_with(
            target_url, data=mock_buffer, headers={'Content-Type': 'application/octet-stream'}, timeout=30
        )
        # 验证是否没有报错
        mock_response.raise_for_status.assert_called_once()

    @patch('requests.Session')
    def test_upload_retry_then_success(self, mock_session_cls, mock_buffer, target_url):
        """测试场景2：第一次失败，第二次成功（验证重试和seek(0)）"""

        # 1. 配置 Mock
        mock_adapter = mock_session_cls.return_value.__enter__.return_value

        # 定义 side_effect：第一次抛出连接错误，第二次返回成功对象
        success_response = MagicMock()
        success_response.status_code = 200

        mock_adapter.put.side_effect = [
            requests.exceptions.ConnectionError("Network down"),  # 第1次：报错
            success_response,  # 第2次：成功
        ]

        # 为了验证 seek(0) 是否被调用，我们需要监控 buffer 对象
        # 使用 MagicMock 包装真实的 BytesIO，或者直接 spy
        with patch.object(mock_buffer, 'seek', wraps=mock_buffer.seek) as mock_seek:
            # 2. 执行函数
            upload_file(url=target_url, buffer=mock_buffer, max_retries=3)

            # 3. 验证
            # 应该调用了2次 put
            assert mock_adapter.put.call_count == 2
            # 应该调用了2次 seek(0)（每次尝试前都会调用）
            assert mock_seek.call_count == 2
            mock_seek.assert_has_calls([call(0), call(0)])

    @patch('requests.Session')
    def test_upload_max_retries_exceeded(self, mock_session_cls, mock_buffer, target_url):
        """测试场景3：重试次数耗尽，抛出异常"""

        # 1. 配置 Mock
        mock_adapter = mock_session_cls.return_value.__enter__.return_value

        # 每次都抛出 Timeout 异常
        mock_adapter.put.side_effect = requests.exceptions.Timeout("Time out")

        # 2. 执行并验证异常
        with pytest.raises(requests.exceptions.Timeout):
            upload_file(url=target_url, buffer=mock_buffer, max_retries=3)

        # 3. 验证调用次数
        # max_retries是3，所以应该尝试了3次
        assert mock_adapter.put.call_count == 3

    @patch('requests.Session')
    def test_upload_500_error_retry(self, mock_session_cls, mock_buffer, target_url):
        """测试场景4：服务器返回500，raise_for_status触发重试"""

        mock_adapter = mock_session_cls.return_value.__enter__.return_value

        # 模拟一个 500 的响应对象
        error_response = MagicMock()
        error_response.status_code = 500
        # 当调用 raise_for_status 时抛出 HTTPError
        error_response.raise_for_status.side_effect = requests.exceptions.HTTPError("500 Server Error")

        # 模拟一个成功的响应对象
        success_response = MagicMock()
        success_response.status_code = 200

        # 流程：第一次500，第二次200
        mock_adapter.put.side_effect = [error_response, success_response]

        upload_file(url=target_url, buffer=mock_buffer, max_retries=3)

        # 验证调用了2次
        assert mock_adapter.put.call_count == 2
