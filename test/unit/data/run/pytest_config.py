import math
import yaml
from swanlab.data.run.main import SwanLabRun, get_run, swanlog
from swanlab.data.run.config import SwanLabConfig, parse, Line, RuntimeInfo
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


def test_parse():
    """
    测试config.parse函数
    """
    # ---------------------------------- omegaConf ----------------------------------
    config_data = {
        "a": 1,
        "b": "mnist",
        "c/d": [1, 2, 3],
        "e/f/h": {"a": 1, "b": {"c": 2}},
    }
    cfg = omegaconf.OmegaConf.create(config_data)
    config = parse(cfg)
    assert yaml.dump(config) == yaml.dump(config_data)
    # ---------------------------------- 包含NaN或者INF ----------------------------------
    config_data = {
        "inf": math.inf,
        "nan": math.nan,
    }
    config = parse(config_data)
    assert config["inf"] == Line.inf
    assert config["nan"] == Line.nan


class TestSwanLabRunConfigOperation:
    """
    单独测试TestSwanLabRunConfig这个类
    """

    def test_basic_operation_object(self):
        """
        测试类的基本操作，增删改
        """
        config = SwanLabConfig()
        # ---------------------------------- 对象风格设置 ----------------------------------
        config.a = 1
        assert config.a == 1
        config.a = 2
        assert config.a == 2
        with pytest.raises(AttributeError):
            config.__a = 1
        # 不存在的属性
        with pytest.raises(AttributeError):
            config.b  # noqa
        # 删除属性
        del config.a
        with pytest.raises(AttributeError):
            del config.a  # 重复删除报错
        with pytest.raises(AttributeError):
            config.a  # noqa

    def test_basic_operation_dict(self):
        """
        测试字典的基本操作，增删改
        """
        config = SwanLabConfig()
        # ---------------------------------- 字典风格设置 ----------------------------------
        config["a"] = 1
        assert config["a"] == 1
        config["a"] = 2
        assert config["a"] == 2
        # 删除属性
        del config["a"]
        with pytest.raises(KeyError):
            del config["a"]
        with pytest.raises(KeyError):
            config["a"]  # noqa
        # 字典风格可以设置，读取，删除私有属性
        config["__a"] = 1
        assert config["__a"] == 1
        del config["__a"]
        with pytest.raises(KeyError):
            del config["__a"]  # 重复删除失败
        with pytest.raises(KeyError):
            config["__a"]  # noqa

    def test_dict_iter(self):
        """
        测试字典风格的迭代
        """
        config = SwanLabConfig()
        ll = ["a", "b", "c", "d"]
        for i in ll:
            config[i] = i
        assert set(config) == {"a", "b", "c", "d"}
        index = 0
        # 返回顺序相同
        for key in config:
            assert key == ll[index]
            index += 1

    def test_dict_len(self):
        """
        测试字典风格的长度
        """
        config = SwanLabConfig()
        assert len(config) == 0
        config["a"] = 1
        assert len(config) == 1
        config["b"] = 2
        assert len(config) == 2
        del config["a"]
        assert len(config) == 1
        del config["b"]
        assert len(config) == 0

    def test_func_operation(self):
        """
        测试内置函数操作
        """
        config = SwanLabConfig()
        # ---------------------------------- get ----------------------------------
        a = config.get("a")
        assert a is None
        a = config.get("a", 1)
        assert a == 1
        config["a"] = 5
        a = config.get("a")
        assert a == 5
        # ---------------------------------- set ----------------------------------
        config.set("b", 1)
        assert config["b"] == 1
        config.set("__b", 1)
        assert config["__b"] == 1
        # ---------------------------------- pop ----------------------------------
        config["c"] = 9
        c = config.pop("c")
        assert c == 9
        assert config.pop("d") is None
        # ---------------------------------- clean ----------------------------------
        config["e"] = 1
        config["g"] = 0
        config.clean()
        assert len(config) == 0
        with pytest.raises(KeyError):
            config["e"]  # noqa
        # ---------------------------------- update ----------------------------------
        config["x"] = 1
        config["y"] = 2
        config["z"] = {"a": 1, "b": 2}
        config.update({"x": 2, "y": 3, "z": {"a": 2, "b": 3}})
        assert config["x"] == 2
        assert config["y"] == 3
        assert config["z"] == {"a": 2, "b": 3}


def test_on_setter():
    """
    测试on_setter函数，在设置属性时触发
    """
    num = 1

    def on_setter(_: RuntimeInfo):
        nonlocal num
        num += 1

    config = SwanLabConfig(on_setter=on_setter)

    # ---------------------------------- 对象、字典风格 ----------------------------------

    # 设置触发
    config.a = 1
    assert num == 2
    del config.a
    assert num == 3
    config["b"] = 1
    assert num == 4
    del config["b"]
    assert num == 5

    # 读取不触发
    config.x = 1
    assert num == 6
    _ = config.x
    assert num == 6
    config["y"] = 1
    assert num == 7
    _ = config["y"]
    assert num == 7

    # ---------------------------------- api ----------------------------------

    # 设置触发
    config.set("c", 1)
    assert num == 8
    config.pop("c")
    assert num == 9
    config.update({"d": {}})
    assert num == 10
    # 深层设置无法触发
    config.d["e"] = 1
    assert num == 10

    # 读取不触发
    config.get("f", 1)
    config["g"] = 1
    assert num == 11
    config.get("g")
    assert num == 11

    # ---------------------------------- clean以后再设置无法触发 ----------------------------------

    config.clean()
    config.h = 1
    assert num == 11


class TestSwanLabRunConfigUseRun:
    """
    测试SwanLabRun的config属性
    """

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
