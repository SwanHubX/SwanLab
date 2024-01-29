from swanlab.utils import get_package_version

def test_get_package_version():
    assert get_package_version() != "unknow"

