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


def test_get_user_setting_path():
    """
    测试获取用户设置文件路径
    """
    assert P.get_user_setting_path() == P.get_host_web() + "/space/~/settings"


def test_get_project_url():
    """
    测试获取项目url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    assert P.get_project_url(username, projname) == P.get_host_web() + "/@" + username + "/" + projname


def test_get_experiment_url():
    """
    测试获取实验url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    expid = nanoid.generate()
    assert (
        P.get_experiment_url(username, projname, expid)
        == P.get_host_web() + "/@" + username + "/" + projname + "/runs/" + expid
    )


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
        测试重复保存，此时会略过
        """
        path = os.path.join(get_save_dir(), ".netrc")
        password = nanoid.generate()
        host = P.get_host_api()
        P.save_key("user", password, host=host)
        change_time = os.path.getmtime(path)
        assert self.get_key(path, host) == password
        P.save_key("user", password, host=host)
        assert self.get_key(path, host) == password
        assert os.path.getmtime(path) == change_time
        time.sleep(0.1)
        # 再次保存，但是账号不同
        new_password = nanoid.generate()
        P.save_key("user2", new_password, host=host)
        assert self.get_key(path, host) == new_password
        assert os.path.getmtime(path) != change_time
        time.sleep(0.1)
        # 再次保存，但是host不同
        new_host = nanoid.generate()
        P.save_key("user", new_password, host=new_host)
        nrc = netrc.netrc(path)
        assert len(nrc.hosts) == 1
        assert nrc.authenticators(new_host) is not None
        time.sleep(0.1)
        # 再次保存，但是密码不同
        new_password = nanoid.generate()
        P.save_key("user", new_password, host=new_host)
        nrc = netrc.netrc(path)
        assert len(nrc.hosts) == 1
        assert nrc.authenticators(new_host)[2] == new_password


class TestIsLogin:
    @staticmethod
    def login():
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
        self.login()
        assert P.is_login()

    def test_no_file(self):
        """
        文件不存在
        """
        assert not P.is_login()

    def test_wrong_host(self):
        """
        host不匹配
        """
        self.login()
        os.environ[SwanLabEnv.API_HOST.value] = nanoid.generate()
        assert not P.is_login()
