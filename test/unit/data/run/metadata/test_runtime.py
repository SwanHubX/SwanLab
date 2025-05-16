import sys

from swanlab import Settings
from swanlab.data.run.metadata.runtime import parse_git_url, get_python_info
import random
import nanoid

from swanlab.package import get_key
from swanlab.swanlab_settings import set_settings


def test_parse_git_url():
    # ssh
    assert parse_git_url("git@github.com:swanhubx/swanlab.git") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("git@localhost:8000/swanhubx/swanlab.git") == "https://localhost:8000/swanhubx/swanlab"
    # https
    assert parse_git_url("https://github.com/swanhubx/swanlab.git") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("https://localhost:8000/swanhubx/swanlab.git") == "https://localhost:8000/swanhubx/swanlab"
    # no .git
    assert parse_git_url("git@github.com:swanhubx/swanlab") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("git@localhost:8000/swanhubx/swanlab") == "https://localhost:8000/swanhubx/swanlab"
    assert parse_git_url("https://github.com/swanhubx/swanlab") == "https://github.com/swanhubx/swanlab"
    assert parse_git_url("https://localhost:8000/swanhubx/swanlab") == "https://localhost:8000/swanhubx/swanlab"

def test_mask_api_key():
    # 1 - 100 之间随机选取数组长度
    length = random.randint(1, 100)
    api_key = get_key()
    # 生成一个长度为 length 的随机字符串数组
    args = [nanoid.generate(size=random.randint(1, 20)) for _ in range(length)]
    sys.argv = args
    # 不含有 api key，不替换
    cmd = get_python_info()["command"]
    assert ("****" not in cmd) is True
    # 存在 api key，替换
    args[random.randint(0, length - 1)] = api_key
    cmd = get_python_info()["command"]
    assert ("****" in cmd) is True
    # 设置为不主动替换隐私信息，即使有 api key 也不替换
    set_settings(Settings(security_mask=False))
    cmd = get_python_info()["command"]
    assert (api_key in cmd) is True
    assert ("****" not in cmd) is True