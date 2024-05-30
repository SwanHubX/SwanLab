from swanlab.package import get_package_latest_version


def test_package_latest_version():
    """
    测试获取swanlab最新版本号
    """
    # 正常情况
    assert isinstance(get_package_latest_version(), str)
    # 超时情况
    assert get_package_latest_version(timeout=1e-3) is None
