from tutils.config import nanoid
import swanlab.package as P
import json
import os

PACKAGE_PATH = os.environ["SWANLAB_PACKAGE_PATH"]

package_data = json.load(open(PACKAGE_PATH))


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


def test_get_host_web():
    """
    测试获取web地址
    """
    assert P.get_host_web() == package_data["host"]["web"]


def test_get_host_api():
    """
    测试获取api地址
    """
    assert P.get_host_api() == package_data["host"]["api"]


def test_get_user_setting_path():
    """
    测试获取用户设置文件路径
    """
    assert P.get_user_setting_path() == P.get_host_web() + "/settings"


def test_get_project_url():
    """
    测试获取项目url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    assert P.get_project_url(username, projname) == P.get_host_web() + "/" + username + "/" + projname


def test_get_experiment_url():
    """
    测试获取实验url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    expid = nanoid.generate()
    assert P.get_experiment_url(
        username, projname,
        expid
    ) == P.get_host_web() + "/" + username + "/" + projname + "/" + expid

# TODO 有关key的测试
