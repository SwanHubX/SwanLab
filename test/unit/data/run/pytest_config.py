from swanlab.data.run.main import SwanLabRun, get_run, SwanLabConfig, swanlog
from tutils import clear, TEMP_PATH
import pytest
import swanlab
import omegaconf


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
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }
        run = SwanLabRun(config=config_data)
        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 4

        assert run.config == config_data

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert run.config["e/f/h"]["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

        assert run.config.a == 1
        assert run.config.b == "mnist"

        assert run.config.get("a") == 1
        assert run.config.get("b") == "mnist"
        assert run.config.get("c/d") == [1, 2, 3]
        assert run.config.get("c/d")[0] == 1
        assert run.config.get("e/f/h") == {"a": 1, "b": {"c": 2}}
        assert run.config.get("e/f/h")["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

    def test_config_null_update(self):
        """
        测试在init之后通过update的方式添加config参数
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        run = SwanLabRun()
        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 0

        run.config.update(config_data)
        assert run.config == config_data
        assert len(run.config) == 4

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert run.config["e/f/h"]["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

        assert run.config.a == 1
        assert run.config.b == "mnist"

        assert run.config.get("a") == 1
        assert run.config.get("b") == "mnist"
        assert run.config.get("c/d") == [1, 2, 3]
        assert run.config.get("c/d")[0] == 1
        assert run.config.get("e/f/h") == {"a": 1, "b": {"c": 2}}
        assert run.config.get("e/f/h")["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

    def test_config_have_update(self):
        """
        测试在init之后通过update的方式添加config参数
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        run = SwanLabRun(config=config_data)
        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 4

        update_data = {
            "a": 2,
            "e/f/h": [4, 5, 6],
            "j": 3,
        }

        run.config.update(update_data)
        assert len(run.config) == 5

        assert run.config["a"] == 2
        assert run.config["e/f/h"] == [4, 5, 6]
        assert run.config["e/f/h"][0] == 4
        assert run.config["j"] == 3

    def test_config_null_set(self):
        """
        测试在没有config初始化时，init之后更新config的参数
        """
        run = SwanLabRun()

        run.config.a = 1
        run.config.b = "mnist"
        run.config["c/d"] = [1, 2, 3]
        run.config.set("e/f/h", {"a": 1, "b": {"c": 2}})

        assert len(run.config) == 4
        assert run.config.a == 1
        assert run.config.b == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}

    def test_config_have_set(self):
        """
        测试在有config初始化时，init之后更新config的参数
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }
        run = SwanLabRun(config=config_data)

        run.config.a = 2
        run.config.set("e/f/h", [4, 5, 6])
        run.config["j"] = 3

        assert len(run.config) == 5
        assert run.config["a"] == 2
        assert run.config["e/f/h"] == [4, 5, 6]
        assert run.config["e/f/h"][0] == 4
        assert run.config["j"] == 3

    def test_config_before_init(self):
        """
        测试在init之前设置的情况
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        swanlab.config.update(config_data)
        run = SwanLabRun()

        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 4

        assert run.config == config_data

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert run.config["e/f/h"]["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

    def test_config_from_omegaconf(self):
        """
        测试config导入omegaconf的情况
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }
        config = omegaconf.OmegaConf.create(config_data)
        run = SwanLabRun(config=config)

        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 4

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == str([1, 2, 3])
        assert run.config["e/f/h"] == str({"a": 1, "b": {"c": 2}})
