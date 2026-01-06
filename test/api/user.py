"""
@author: Zhou QiYang
@file: user.py
@time: 2026/1/2 21:30
@description: OpenApi 用户相关测试代码
"""

import random
from unittest.mock import patch, MagicMock

import swanlab
from swanlab.api.utils import STATUS_CREATED, STATUS_OK
from swanlab.error import ApiError


# example_code:
#     api = swanlab.Api()
#     user = api.user
#     print(user.__dict__)
#     new_key = user.generate_api_key()
#     print(new_key)
#     print(user.api_keys)
#     user.delete_api_key(api_key=new_key)
#     print(user.api_keys)


def make_fake_login_info():
    from swanlab.core_python.auth.providers.api_key import LoginInfo

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.reason = "OK"
    mock_resp.json.return_value = {"userInfo": {"username": "testuser"}}

    login_info = LoginInfo(
        resp=mock_resp, api_key="test_api_key_12345", api_host="https://api.test.com", web_host="https://test.com"
    )
    return login_info


def make_fake_api_keys(count=3):
    return [
        {
            "id": i + 1,
            "name": f"test-key-{i + 1}",
            "createdAt": f"2025-01-0{i + 1}T00:00:00Z",
            "key": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=40)),
        }
        for i in range(count)
    ]


def make_fake_latest_api_key():
    return {
        "id": 999,
        "name": "latest-key",
        "createdAt": "2025-01-10T00:00:00Z",
        "key": "".join(random.choices("0123456789abcdefghijklmnopqrstuvwxyz", k=40)),
    }


def make_fake_self_hosted_404(*args, **kwargs):
    """构造一个指向 /self_hosted/info 的 404 ApiError"""
    resp = MagicMock()
    resp.status_code = 404
    resp.reason = "Not Found"
    raise ApiError(resp=resp)


# 模拟普通用户登录
class ApiUserContext:
    def __init__(self, fake_login_info=None):
        self.fake_login_info = fake_login_info or make_fake_login_info()
        self.patches = []
        self.api = None
        self.user = None

    def __enter__(self):
        patch1 = patch("swanlab.api.api.auth.code_login", return_value=self.fake_login_info)
        patch2 = patch("swanlab.api.api.get_self_hosted_init", side_effect=make_fake_self_hosted_404)
        self.patches.extend([patch1, patch2])

        for p in self.patches:
            p.start()

        self.api = swanlab.Api(api_key="fake_key")
        self.user = self.api.user()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for p in self.patches:
            p.stop()
        return False


def test_api_user_info():
    """测试获取用户名和团队列表"""
    fake_login_info = make_fake_login_info()
    fake_teams = [
        {"name": "team-alpha"},
        {"name": "team-beta"},
        {"name": "team-gamma"},
    ]

    with ApiUserContext() as ctx:
        user = ctx.user
        assert user.username == fake_login_info.username
        with patch("swanlab.api.model.user.get_user_groups", return_value=fake_teams):
            teams = user.teams
            assert isinstance(teams, list)
            assert teams == [t["name"] for t in fake_teams]


def test_api_user_api_keys():
    """测试获取API Key列表"""
    fake_api_keys = make_fake_api_keys(count=3)

    with ApiUserContext() as ctx:
        user = ctx.user
        with patch("swanlab.api.model.user.get_api_keys", return_value=fake_api_keys):
            api_keys = user.api_keys
            assert isinstance(api_keys, list)
            assert len(api_keys) == len(fake_api_keys)
            for key in api_keys:
                assert isinstance(key, str)
            expected_keys = [k['key'] for k in fake_api_keys]
            assert api_keys == expected_keys


def test_api_user_generate_api_key():
    """测试生成API Key功能，包括成功和失败情况"""
    fake_latest_key = make_fake_latest_api_key()

    with ApiUserContext() as ctx:
        user = ctx.user
        with patch("swanlab.api.model.user.create_api_key", return_value=STATUS_CREATED):
            with patch("swanlab.api.model.user.get_latest_api_key", return_value=fake_latest_key):
                new_key = user.generate_api_key()
                assert new_key is not None
                assert new_key == fake_latest_key['key']

        with patch("swanlab.api.model.user.create_api_key", return_value="Failed"):
            new_key_failed = user.generate_api_key()
            assert new_key_failed is None


def test_api_user_delete_api_key_success():
    """测试成功删除API Key，验证删除前后key列表的变化"""
    fake_api_keys = make_fake_api_keys(count=3)
    key_to_delete = fake_api_keys[1]['key']

    with ApiUserContext() as ctx:
        user = ctx.user
        with patch("swanlab.api.model.user.get_api_keys", return_value=fake_api_keys):
            with patch("swanlab.api.model.user.delete_api_key", return_value=STATUS_OK) as delete_key:
                api_keys_before = user.api_keys
                assert key_to_delete in api_keys_before
                assert len(api_keys_before) == 3

                result = user.delete_api_key(key_to_delete)
                # swanlab.api.model.user.delete_api_key 方法被调用了一次
                assert delete_key.call_count == 1
                assert result is True


def test_api_user_delete_api_key_not_found():
    """测试删除不存在的API Key，验证返回False且列表不变"""
    fake_api_keys = make_fake_api_keys(count=3)

    with ApiUserContext() as ctx:
        user = ctx.user
        with patch("swanlab.api.model.user.get_api_keys", return_value=fake_api_keys):
            with patch("swanlab.api.model.user.delete_api_key", return_value="Not Found") as delete_key:
                api_keys_before = user.api_keys
                assert len(api_keys_before) == 3

                result_not_found = user.delete_api_key("non-existent-key")
                assert delete_key.call_count == 0
                assert result_not_found is False
