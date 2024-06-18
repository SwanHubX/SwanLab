from tutils.config import nanoid
import json
from swanlab.package import (
    get_package_version,
    get_host_web,
    get_user_setting_path,
    get_host_api,
    get_project_url,
    get_experiment_url
)
import os
from swanlab.package import get_package_latest_version

PACKAGE_PATH = os.environ["SWANLAB_PACKAGE_PATH"]

package_data = json.load(open(PACKAGE_PATH))


def test_package_latest_version():
    """
    测试获取swanlab最新版本号
    """
    # 正常情况
    assert isinstance(get_package_latest_version(), str)
    # 超时情况
    assert get_package_latest_version(timeout=1e-3) is None


def mock_get_package_version():
    return get_package_version(PACKAGE_PATH)


def mock_get_host_web():
    return get_host_web(PACKAGE_PATH)


def mock_get_host_api():
    return get_host_api(PACKAGE_PATH)


def mock_get_user_setting_path():
    return get_user_setting_path(PACKAGE_PATH)


def mock_get_project_url(username: str, projname: str):
    return get_project_url(username, projname, PACKAGE_PATH)


def mock_get_experiment_url(username: str, projname: str, expid: str):
    return get_experiment_url(username, projname, expid, PACKAGE_PATH)


# ---------------------------------- 简单测试一下 ----------------------------------


def test_get_package_version():
    """
    测试获取版本号
    """
    assert mock_get_package_version() == package_data["version"]


def test_get_host_web():
    """
    测试获取web地址
    """
    assert mock_get_host_web() == package_data["host"]["web"]


def test_get_host_api():
    """
    测试获取api地址
    """
    assert mock_get_host_api() == package_data["host"]["api"]


def test_get_user_setting_path():
    """
    测试获取用户设置文件路径
    """
    assert mock_get_user_setting_path() == mock_get_host_web() + "/settings"


def test_get_project_url():
    """
    测试获取项目url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    assert mock_get_project_url(username, projname) == mock_get_host_web() + "/" + username + "/" + projname


def test_get_experiment_url():
    """
    测试获取实验url
    """
    username = nanoid.generate()
    projname = nanoid.generate()
    expid = nanoid.generate()
    assert (mock_get_experiment_url(username, projname, expid)
            == mock_get_host_web() + "/" + username + "/" + projname + "/" + expid)
