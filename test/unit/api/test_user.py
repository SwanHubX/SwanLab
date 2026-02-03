from unittest.mock import patch, MagicMock

import pytest

from swanlab.api.user import User
from swanlab.api.users import Users
from swanlab.core_python import Client
from swanlab.error import ApiError
from utils import create_user_data


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


@pytest.mark.parametrize(
    "kwargs",
    [
        {},  # 私有化未启动
        {'start': True, 'enabled': False},  # 启动，但未启用
        {'start': True, 'expired': True},  # licence过期
        {'start': True},  # 启动，正常，但不是root
        {'username': 'test_other_user'},
    ],
)
def test_create_permission(kwargs):
    """测试尝试创建用户是否会被拦截"""
    user = create_user(kwargs.get('username', None))
    with SelfHosted(**kwargs):
        with pytest.raises(ValueError):
            user.create(username='test_user', password='123456aa')


@pytest.mark.parametrize(
    ("username", "password"),
    [
        ('user@name', 'password123'),  # 无效的用户名
        ('test_user', 'short'),  # 无效密码（密码长度小于8）
        ('test_user', '12345678'),  # 无效密码（全是数字）
        ('test_user', 'ABCDEFGH'),  # 有效密码（全是字母）
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
        with pytest.raises(ValueError):
            assert other_user.generate_api_key() is None
        with pytest.raises(ValueError):
            assert other_user.delete_api_key(api_key='test_api_key') == False


def test_users():
    """测试能否分页获取所有用户"""
    with patch('swanlab.api.users.get_users') as mock_get_users:
        total = 80
        page_size = 20

        def side_effect(*args, **kwargs):
            return create_user_data(page=kwargs.get("page", 1), total=total)

        mock_get_users.side_effect = side_effect
        client = MagicMock(spec=Client)
        users = Users(client, login_user="test_user")

        user_list = list(users)
        assert len(user_list) == total
        for i, user in enumerate(user_list):
            assert user.username == f'user_{i}'

        assert mock_get_users.call_count == (total + page_size - 1) // page_size
