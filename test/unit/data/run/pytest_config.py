from swanlab.data.run.main import SwanLabRun, get_run, SwanLabConfig, swanlog
from tutils import clear
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
    def test_config_have(self):
        """
        初始化时有config参数，测试三种获取数据的方式
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }
        run = SwanLabRun(run_config=config_data)
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
        测试init为空，之后通过update的方式添加config参数
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

    def test_config_have_update(self):
        """
        测试init不为空，之后通过update的方式添加config参数
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        run = SwanLabRun(run_config=config_data)
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
        run = SwanLabRun(run_config=config_data)

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

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["c/d"][0] == 1
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert run.config["e/f/h"]["a"] == 1
        assert run.config["e/f/h"]["b"]["c"] == 2

    def test_config_clean(self):
        """
        测试在finish之后config是否置空
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        run = SwanLabRun(run_config=config_data)
        run.finish()

        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 0

    def test_config_get_config(self):
        """
        测试get_config
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        assert isinstance(swanlab.get_config(), SwanLabConfig)
        assert len(swanlab.get_config()) == 0

        run = SwanLabRun(run_config=config_data)

        assert isinstance(swanlab.get_config(), SwanLabConfig)
        assert len(swanlab.get_config()) == 4

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
        cfg = omegaconf.OmegaConf.create(config_data)
        run = SwanLabRun(run_config=cfg)

        assert isinstance(run.config, SwanLabConfig)
        assert len(run.config) == 4

        assert run.config["a"] == 1
        assert run.config["b"] == "mnist"
        assert run.config["c/d"] == [1, 2, 3]
        assert run.config["e/f/h"] == {"a": 1, "b": {"c": 2}}

    def test_not_json_serializable(self):
        """
        测试不可json化的数据
        """
        import math, json

        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
            "test_nan": math.nan,
            "test_inf": math.inf,
        }

        run = SwanLabRun(run_config=config_data)

        assert run.config.test_nan == "nan"
        assert run.config.test_inf == "inf"

        json_data = json.dumps(dict(run.config))

    # def test_insert_class(self):
    #     """
    #     测试插入类
    #     """
    #     from collections.abc import MutableMapping

    #     class Test(MutableMapping):
    #         def __init__(self, a, b):
    #             self.data = {"a": a, "b": b}

    #         def __setitem__(self, __key, __value):
    #             self.data[__key] = __value

    #         def __delitem__(self, __key):
    #             del self.data[__key]

    #         def __getitem__(self, __key):
    #             return self.data.get(__key, None)

    #         def __len__(self):
    #             return len(self.data)

    #         def __iter__(self):
    #             return iter(self.data)

    #     config_data = {"a": 1, "b": "mnist", "c/d": [1, 2, 3], "e/f/h": {"a": 1, "b": {"c": 2}}, "test": Test(1, 2)}

    #     run = SwanLabRun(run_config=config_data)
    #     assert run.config.test.a == 1
    #     assert run.config.test.b == 2
