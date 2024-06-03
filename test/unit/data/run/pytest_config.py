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
    def test_config_normal_haverun(self):
        """
        初始化时有config参数，测试三种获取数据的方式，且使用的run对象
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

        run.config.save()

    def test_config_finish_haverun(self):
        """
        测试在run.finish()之后config是否置空
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

        run.config.save()

    def test_config_normal(self):
        """
        初始化时有config参数，测试三种获取数据的方式，直接用全局的config对象
        """
        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }
        config = SwanLabConfig(config=config_data)

        assert isinstance(config, SwanLabConfig)
        assert len(config) == 4

        assert config == config_data

        assert config["a"] == 1
        assert config["b"] == "mnist"
        assert config["c/d"] == [1, 2, 3]
        assert config["c/d"][0] == 1
        assert config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert config["e/f/h"]["a"] == 1
        assert config["e/f/h"]["b"]["c"] == 2

        assert config.a == 1
        assert config.b == "mnist"

        assert config.get("a") == 1
        assert config.get("b") == "mnist"
        assert config.get("c/d") == [1, 2, 3]
        assert config.get("c/d")[0] == 1
        assert config.get("e/f/h") == {"a": 1, "b": {"c": 2}}
        assert config.get("e/f/h")["a"] == 1
        assert config["e/f/h"]["b"]["c"] == 2

        config.save()
        config.clean()

    def test_config_update(self):
        """
        测试config初始为空，之后通过update的方式添加config参数
        """

        config_data = {
            "a": 1,
            "b": "mnist",
            "c/d": [1, 2, 3],
            "e/f/h": {"a": 1, "b": {"c": 2}},
        }

        update_data = {
            "a": 2,
            "e/f/h": [4, 5, 6],
            "j": 3,
        }

        config = SwanLabConfig()
        assert len(config) == 0

        # 第一次更新
        config.update(config_data)
        assert config == config_data
        assert len(config) == 4

        assert config["a"] == 1
        assert config["b"] == "mnist"
        assert config["c/d"] == [1, 2, 3]
        assert config["c/d"][0] == 1
        assert config["e/f/h"] == {"a": 1, "b": {"c": 2}}
        assert config["e/f/h"]["a"] == 1
        assert config["e/f/h"]["b"]["c"] == 2

        # 第二次更新
        config.update(update_data)
        assert len(config) == 5

        assert config["a"] == 2
        assert config["e/f/h"] == [4, 5, 6]
        assert config["e/f/h"][0] == 4
        assert config["j"] == 3

        config.save()
        config.clean()

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

        config = SwanLabConfig(config=config_data)

        assert isinstance(swanlab.get_config(), SwanLabConfig)
        assert len(swanlab.get_config()) == 4

        config.save()
        config.clean()

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
        config = SwanLabConfig(config=config_data)

        assert isinstance(config, SwanLabConfig)
        assert len(config) == 4

        assert config["a"] == 1
        assert config["b"] == "mnist"
        assert config["c/d"] == [1, 2, 3]
        assert config["e/f/h"] == {"a": 1, "b": {"c": 2}}

        config.save()
        config.clean()

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

        config = SwanLabConfig(config=config_data)

        json_data = json.dumps(dict(config))

        config.save()
        config.clean()

    def test_insert_class(self):
        """
        测试插入类
        """
        from collections.abc import MutableMapping

        class Test(MutableMapping):
            def __init__(self, a, b):
                self.data = {"a": a, "b": b}

            def __setitem__(self, __key, __value):
                self.data[__key] = __value

            def __delitem__(self, __key):
                del self.data[__key]

            def __getitem__(self, __key):
                return self.data.get(__key, None)

            def __len__(self):
                return len(self.data)

            def __iter__(self):
                return iter(self.data)

        config_data = {"a": 1, "b": "mnist", "c/d": [1, 2, 3], "e/f/h": {"a": 1, "b": {"c": 2}}, "test": Test(1, 2)}

        config = SwanLabConfig(config=config_data)

        assert config.test.data["a"] == 1
        assert config.test.data["b"] == 2

        config.save()
        config.clean()
