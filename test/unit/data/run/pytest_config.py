from swanlab.data.run.main import SwanLabRun, get_run, SwanLabRunState, swanlog
from tutils import clear, TEMP_PATH
import pytest


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    if get_run() is not None:
        get_run().finish()
    swanlog.disable_log()
    yield
    clear()
    swanlog.enable_log()


class TestSwanLabRunConfig:
    def test_config_ok(self):
        """
        测试解析一个正常的数字
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": 2},
        }
        run = SwanLabRun(config=config_data)
        assert run.config == config_data

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": 2}
        assert run.config["e/f/h"]["a"] == 1

        assert run.config.a == 1
        assert run.config.b == "mnist"

        assert run.config.get("a") == 1
        assert run.config.get("b") == "mnist"
        assert run.config.get("c/d") == [1, 2, 3]
        assert run.config.get("c/d")[0] == 1
        assert run.config.get("e/f/h") == {"a": 1, "b": 2}
        assert run.config.get("e/f/h")["a"] == 1

    def test_config_update(self):
        """
        测试在init之后通过update的方式添加config参数
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": 2},
        }

        run = SwanLabRun()
        # assert len(run.config) == 0

        run.config.update(config_data)
        assert run.config == config_data

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": 2}
        assert run.config["e/f/h"]["a"] == 1

        assert run.config.a == 1
        assert run.config.b == "mnist"

        assert run.config.get("a") == 1
        assert run.config.get("b") == "mnist"
        assert run.config.get("c/d") == [1, 2, 3]
        assert run.config.get("c/d")[0] == 1
        assert run.config.get("e/f/h") == {"a": 1, "b": 2}
        assert run.config.get("e/f/h")["a"] == 1
