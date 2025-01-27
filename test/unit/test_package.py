import json
import netrc
import os
import time

import nanoid
import pytest

import swanlab.package as P
from swanlab.env import SwanLabEnv, get_save_dir

_ = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
package_data = json.load(open(os.path.join(_, "swanlab", "package.json")))


def test_package_latest_version():
    """
    测试获取swanlab最新版本号
    """
    # 正常情况
    assert isinstance(P.get_package_latest_version(), str)
    # 超时情况
    assert P.get_package_latest_version(timeout=1e-3) is None


def test_get_package_version():
    """
    测试获取版本号
    """
    assert P.get_package_version() == package_data["version"]


def test_get_host_web_env():
    """
    通过环境变量指定web地址
    """
    os.environ[SwanLabEnv.WEB_HOST.value] = nanoid.generate()
    assert P.get_host_web() == os.environ[SwanLabEnv.WEB_HOST.value]


def test_get_host_api_env():
    """
    通过环境变量指定api地址
    """
    os.environ[SwanLabEnv.API_HOST.value] = nanoid.generate()
    assert P.get_host_api() == os.environ[SwanLabEnv.API_HOST.value]


def test_fmt_web_host():
    """
    测试格式化web地址
    """
    assert P.fmt_web_host() == P.get_host_web().rstrip("/")
    assert P.fmt_web_host("https://abc.cn/") == "https://abc.cn"


# ---------------------------------- 登录部分 ----------------------------------


class TestGetKey:
    @staticmethod
    def remove_env_key():
        if SwanLabEnv.API_KEY.value in os.environ:
            del os.environ[SwanLabEnv.API_KEY.value]

    def test_ok(self):
        """
        获取key成功
        """
        self.remove_env_key()
        # 首先需要登录
        file = os.path.join(get_save_dir(), ".netrc")
        with open(file, "w"):
            pass
        nrc = netrc.netrc(file)
        key = nanoid.generate()
        nrc.hosts[P.get_host_api()] = ("user", "", key)
        with open(file, "w") as f:
            f.write(nrc.__repr__())
        assert P.get_key() == key

    def test_no_file(self):
        """
        文件不存在
        """
        self.remove_env_key()
        from swanlab.error import KeyFileError

        with pytest.raises(KeyFileError) as e:
            P.get_key()
        assert str(e.value) == "The file does not exist"

    def test_no_host(self):
        from swanlab.error import KeyFileError

        self.test_ok()  # 此时删除了环境变量
        host = nanoid.generate()
        os.environ[SwanLabEnv.API_HOST.value] = host
        assert P.get_host_api() == host
        with pytest.raises(KeyFileError) as e:
            P.get_key()
        assert str(e.value) == f"The host {host} does not exist"

    def test_use_env(self):
        """
        使用环境变量，优先级高于本地文件
        """
        self.test_ok()
        key = nanoid.generate()
        os.environ[SwanLabEnv.API_KEY.value] = key
        assert P.get_key() == key


class TestSaveKey:

    @staticmethod
    def get_key(path, host):
        nrc = netrc.netrc(path)
        info = nrc.authenticators(host)
        return info[2]

    @staticmethod
    def loose_compare(a, b, equal: bool = True):
        """
        不同文件系统的时间精度不同，因此需要 sleep 一段时间进行比较
        不精确到微妙 只精确到毫秒
        """
        assert (abs(a - b) < 1e-5) is equal

    def test_ok(self):
        """
        保存key成功
        """
        path = os.path.join(get_save_dir(), ".netrc")
        password = nanoid.generate()
        host = P.get_host_api()
        P.save_key("user", password, host=host)
        assert self.get_key(path, host) == password
        # 在保存一次，保证只存在一个host
        new_host = nanoid.generate()
        P.save_key("user", password, host=new_host)
        nrc = netrc.netrc(path)
        assert len(nrc.hosts) == 1
        assert nrc.authenticators(new_host) is not None

    def test_duplicate(self):
        """
        测试重复保存，此时会略过保存，因此不会改变文件的修改时间
        """

        path = os.path.join(get_save_dir(), ".netrc")

        def duplicate_save(p: str, h: str, user: str = "user", equal: bool = True):
            c = os.path.getmtime(path)
            time.sleep(0.1)
            P.save_key(user, p, host=h)
            assert self.get_key(path, h) == p
            return self.loose_compare(os.path.getmtime(path), c, equal)

        password = nanoid.generate()
        host = P.get_host_api()
        P.save_key("user", password, host=host)
        assert self.get_key(path, host) == password

        # 重复保存
        duplicate_save(password, host, equal=True)
        # 再次保存，但是账号不同
        new_password = nanoid.generate()
        duplicate_save(new_password, host, user="user2", equal=False)
        # 再次保存，但是host不同
        new_host = nanoid.generate()
        duplicate_save(new_password, new_host, equal=False)
        nrc = netrc.netrc(path)
        assert len(nrc.hosts) == 1
        assert nrc.authenticators(new_host) is not None
        # 再次保存，但是密码又不同了
        new_password = nanoid.generate()
        duplicate_save(new_password, new_host, equal=False)
        nrc = netrc.netrc(path)
        assert len(nrc.hosts) == 1
        assert nrc.authenticators(new_host)[2] == new_password
        # 再次保存，但是密码又相同了
        duplicate_save(new_password, new_host, equal=True)


class TestHasApiKey:
    @staticmethod
    def save_api_key():
        path = os.path.join(get_save_dir(), ".netrc")
        with open(path, "w"):
            pass
        nrc = netrc.netrc(path)
        key = nanoid.generate()
        nrc.hosts[P.get_host_api()] = ("user", "", key)
        with open(path, "w") as f:
            f.write(nrc.__repr__())

    def test_ok(self):
        """
        已经登录
        """
        self.save_api_key()
        assert P.has_api_key()

    def test_no_file(self):
        """
        文件不存在
        """
        assert not P.has_api_key()

    def test_wrong_host(self):
        """
        host不匹配
        """
        self.save_api_key()
        os.environ[SwanLabEnv.API_HOST.value] = nanoid.generate()
        assert not P.has_api_key()
