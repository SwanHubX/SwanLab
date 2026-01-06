"""
@author: Zhou QiYang
@file: self_hosted.py
@time: 2026/1/6 15:28
@description: OpenApi 私有化相关功能测试
"""

from unittest.mock import MagicMock, patch

import pytest

import swanlab
from swanlab.api.model.user import ApiUser, SuperUser
from swanlab.api.utils import STATUS_CREATED


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


# 构造一个 root 商业版私有化部署信息，用于生成 SuperUser
def make_fake_self_hosted_root_info(enabled=True, expired=False, root=True, plan="commercial", seats=99) -> dict:
    return {
        "enabled": enabled,
        "expired": expired,
        "root": root,
        "plan": plan,
        "seats": seats,
    }


def _user_with_self_hosted(self_hosted_info: dict):
    """
    内部工具函数：根据给定 self_hosted_info 构造一个 Api 实例并返回 user 对象，
    不对结果类型做任何断言，方便不同测试场景复用。
    """
    fake_login_info = make_fake_login_info()
    with patch("swanlab.api.api.auth.code_login", return_value=fake_login_info), patch(
        "swanlab.api.api.get_self_hosted_init", return_value=self_hosted_info
    ):
        api = swanlab.Api(api_key="fake_key")
        # 这里传入 username，模拟 self-hosted 环境下按用户名获取用户
        return api.user(username="admin")


@pytest.mark.parametrize(
    "enabled, expired, root, plan, expected_type, expect_exception",
    [
        # 唯一应该返回 SuperUser 的组合
        (True, False, True, "commercial", SuperUser, None),
        # 1) 未启用 -> RuntimeError
        (False, False, True, "commercial", None, RuntimeError),
        # 2) 已过期 -> RuntimeError
        (True, True, True, "commercial", None, RuntimeError),
        # 3) free 版本 -> ApiUser
        (True, False, True, "free", ApiUser, None),
        # 4) commercial 但非 root -> ApiUser
        (True, False, False, "commercial", ApiUser, None),
        # 5) 其他 plan（预留）-> None
        (True, False, True, "education", type(None), None),
    ],
)
def test_super_user_created(enabled, expired, root, plan, expected_type, expect_exception):
    """
    仅当 enabled=True, expired=False, root=True, plan='commercial' 时，
    Api.user 返回 SuperUser；其它组合要么抛出异常，要么返回 ApiUser/None。
    """
    info = make_fake_self_hosted_root_info(enabled=enabled, expired=expired, root=root, plan=plan)

    if expect_exception is not None:
        with pytest.raises(expect_exception):
            _user_with_self_hosted(info)
    else:
        user = _user_with_self_hosted(info)
        assert isinstance(user, expected_type)


@pytest.mark.parametrize(
    "case, username, password, create_ret, expect_exception, expect_result, expect_called",
    [
        # 正常情况：用户名、密码都合法，底层 create_user 返回 STATUS_CREATED
        ("ok", "new_user", "Password123", "STATUS_CREATED", None, True, 1),
        # 用户名非法：包含空格或特殊字符 -> 立刻 ValueError，不调用 create_user
        ("bad_username", "bad user", "Password123", None, ValueError, None, 0),
        # 密码非法：太短或不包含数字/字母 -> 立刻 ValueError，不调用 create_user
        ("bad_password", "good_user", "short", None, ValueError, None, 0),
        # 底层 create_user 未返回 STATUS_CREATED -> 触发 raise False -> TypeError
        ("backend_failed", "new_user", "Password123", "Failed", TypeError, None, 1),
    ],
)
def test_super_user_create_parametrized(
        case, username, password, create_ret, expect_exception, expect_result, expect_called
):
    """
    在 self-hosted root + commercial 环境下，使用参数化方式测试 SuperUser.create：
    - 覆盖用户名/密码校验不通过抛出的 ValueError
    - 覆盖底层 create_user 成功/失败时的行为
    """
    info = make_fake_self_hosted_root_info(enabled=True, expired=False, root=True, plan="commercial")

    fake_login_info = make_fake_login_info()
    with patch("swanlab.api.api.auth.code_login", return_value=fake_login_info), patch(
            "swanlab.api.api.get_self_hosted_init", return_value=info
    ):
        api = swanlab.Api(api_key="fake_key")
        super_user = api.user(username="admin")
        assert isinstance(super_user, SuperUser)

        # 无论是否会真正调用 create_user，一律 patch 掉，便于统计调用次数
        with patch("swanlab.api.model.user.create_user") as mock_create:
            if create_ret is not None:
                # 根据用例配置底层返回值：STATUS_CREATED 或 Failed
                mock_create.return_value = STATUS_CREATED if create_ret == "STATUS_CREATED" else create_ret

            if expect_exception is not None:
                with pytest.raises(expect_exception):
                    super_user.create(username=username, password=password)
                # 断言调用次数
                assert mock_create.call_count == expect_called
            else:
                result = super_user.create(username=username, password=password)
                assert result is expect_result
                assert mock_create.call_count == expect_called
