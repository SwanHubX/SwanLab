from unittest.mock import MagicMock

import pytest

from swanlab.api.user import User
from swanlab.core_python import Client
from swanlab.core_python.api.type import IdentityType


def create_user(username=None, identity: IdentityType = "user"):
    """创建用户对象的辅助函数"""
    return User(MagicMock(spec=Client), login_user="test_user", username=username, identity=identity)


def test_create_permission():
    """测试普通用户尝试创建用户是否会被拦截"""
    user = create_user(identity="user")
    assert user.create(username='test_user', password='123456aa') == False


@pytest.mark.parametrize(
    ("username", "password"),
    [
        ('user@name', 'password123'),
        ('test_user', 'short'),
        ('test_user', '12345678'),
        ('test_user', 'ABCDEFGH'),
    ],
)
def test_check_create_info(username, password):
    """测试无效的用户名或密码"""
    root_user = create_user(identity="root")
    with pytest.raises(ValueError):
        root_user.create(username, password)


def test_other_user():
    """测试是否对未开发的功能进行拦截"""
    other_user = create_user(identity="root", username="other_user")
    assert other_user.generate_api_key() is None
    assert other_user.delete_api_key(api_key='test_api_key') == False
