from unittest.mock import patch, MagicMock

import pytest

from swanlab.api.user import User
from swanlab.core_python import Client
from swanlab.error import ApiError


def create_user(username=None):
    """创建用户对象的辅助函数"""
    return User(MagicMock(spec=Client), login_user="test_user", username=username)


class SelfHosted:
    """私有化部署的上下文管理器"""

    def __init__(self, start=False, **kwargs):
        self.start = start  # 是否启动私有化部署
        self.enabled = kwargs.get('enabled', True)
        self.expired = kwargs.get('expired', False)
        self.root = kwargs.get('root', False)
        self.plan = kwargs.get('plan', 'free')
        self.seats = kwargs.get('seats', 99)

    def __enter__(self):
        self.mock_get_metrics = patch('swanlab.api.utils.get_self_hosted_init').start()
        if not self.start:
            self.mock_get_metrics.side_effect = ApiError()
        else:
            self.mock_get_metrics.return_value = {
                'enabled': self.enabled,  # 是否成功部署
                'expired': self.expired,  # licence是否过期
                'root': self.root,  # 是否为根用户
                'plan': self.plan,  # 私有化版本（免费、商业）
                'seats': self.seats,  # 余剩席位
            }
        return self

    def __exit__(self, *args):
        patch.stopall()


def test_create_permission():
    """测试普通用户尝试创建用户是否会被拦截"""
    user = create_user()
    with SelfHosted():
        with pytest.raises(ValueError):
            user.create(username='test_user', password='123456aa')
    with SelfHosted(start=True, enabled=False):
        assert user.create(username='test_user', password='123456aa') is None
    with SelfHosted(start=True, expired=True):
        assert user.create(username='test_user', password='123456aa') is None
    with SelfHosted(start=True):
        assert user.create(username='test_user', password='123456aa') is None


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
    root_user = create_user()
    with SelfHosted(start=True, root=True):
        with pytest.raises(ValueError):
            root_user.create(username, password)


def test_other_user():
    """测试是否对未开发的功能进行拦截"""
    other_user = create_user(username="other_user")
    with SelfHosted(start=True, root=True):
        assert other_user.generate_api_key() is None
        assert other_user.delete_api_key(api_key='test_api_key') == False
