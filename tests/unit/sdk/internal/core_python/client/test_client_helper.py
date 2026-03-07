"""
@author: cunyue
@file: test_client_helper.py
@time: 2026/3/7 21:55
@description: 测试 SwanLab 运行时客户端辅助函数
"""

from unittest.mock import MagicMock

import requests

from swanlab.sdk.internal.core_python.client.helper import decode_error_response, decode_response


def test_decode_response_json():
    """测试正常解码 JSON 响应"""
    mock_resp = MagicMock(spec=requests.Response)
    mock_resp.json.return_value = {"message": "success"}

    assert decode_response(mock_resp) == {"message": "success"}


def test_decode_response_fallback_text():
    """测试非 JSON 格式响应时，退化为返回 text"""
    mock_resp = MagicMock(spec=requests.Response)
    # 模拟 json() 解析抛出异常
    mock_resp.json.side_effect = requests.JSONDecodeError("Expecting value", "", 0)
    mock_resp.text = "<html>502 Bad Gateway</html>"

    assert decode_response(mock_resp) == "<html>502 Bad Gateway</html>"


def test_decode_error_response_empty():
    """测试空响应体，应直接返回 None"""
    mock_resp = MagicMock(spec=requests.Response)
    mock_resp.text = "   "

    assert decode_error_response(mock_resp) is None


def test_decode_error_response_valid_json():
    """测试正常的带 code 和 message 的 JSON 错误响应"""
    mock_resp = MagicMock(spec=requests.Response)
    mock_resp.text = '{"code": 4001, "message": "Invalid Params"}'
    mock_resp.json.return_value = {"code": 4001, "message": "Invalid Params"}
    mock_resp.status_code = 400
    mock_resp.reason = "Bad Request"

    result = decode_error_response(mock_resp)
    assert result == ("4001", "Invalid Params")


def test_decode_error_response_missing_fields():
    """测试 JSON 格式正确，但缺失 code/message 字段时，使用状态码 fallback"""
    mock_resp = MagicMock(spec=requests.Response)
    mock_resp.text = '{"error": "something went wrong"}'
    mock_resp.json.return_value = {"error": "something went wrong"}
    mock_resp.status_code = 404
    mock_resp.reason = "Not Found"

    result = decode_error_response(mock_resp)
    assert result == ("404", "Not Found")


def test_decode_error_response_invalid_json():
    """测试非 JSON 错误响应，装饰器应捕获异常并返回 None"""
    mock_resp = MagicMock(spec=requests.Response)
    mock_resp.text = "Internal Server Error"
    mock_resp.json.side_effect = requests.JSONDecodeError("Expecting value", "", 0)

    assert decode_error_response(mock_resp) is None
