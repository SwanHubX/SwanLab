import math
import yaml
from swanlab.data.run.main import SwanLabRun, get_run, swanlog, get_config
from swanlab.data.run.config import SwanLabConfig, parse, Line, RuntimeInfo, MutableMapping
import pytest
import omegaconf
from dataclasses import dataclass
import argparse


@pytest.fixture(scope="function", autouse=True)
def setup_function():
    """
    在当前测试文件下的每个测试函数执行前后执行
    """
    if get_run() is not None:
        get_run().finish()
    swanlog.disable_log()
    yield
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

    # ---------------------------------- 自定义继承自MutableMapping的类 ----------------------------------

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
    config = parse(config_data)
    assert config["test"]["a"] == 1
    assert config["test"]["b"] == 2
    # ---------------------------------- 包含NaN或者INF的dict对象 ----------------------------------
    config_data = {
        "inf": math.inf,
        "nan": math.nan,
    }
    config = parse(config_data)
    assert config["inf"] == Line.inf
    assert config["nan"] == Line.nan

    # ---------------------------------- dataclass support ----------------------------------
    @dataclass
    class MyData:
        a: int
        b: float

    config_data = MyData(10, 20.0)
    config = parse(config_data)
    assert config["a"] == 10
    assert config["b"] == 20.0
    # ---------------------------------- argparse.Namespace ----------------------------------
    config_data = argparse.Namespace(a=1, b="mnist", c=[1, 2, 3], d={"a": 1, "b": {"c": 2}})
    config = parse(config_data)
    assert yaml.dump(config) == yaml.dump(vars(config_data))


def test_parse_base_class():
    """
    继承自基础的类不能绕过parse函数
    """

    class StrChild(str):
        pass

    value = parse(StrChild("abc"))
    assert value == "abc"
    assert yaml.safe_dump({"value": value})

    class IntChild(int):
        pass

    value = parse(IntChild(1))
    assert value == 1
    assert value != "1"
    assert yaml.safe_dump({"value": value})


class TestSwanLabConfigOperation:
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
        # int访问，设置
        config[1] = 1  # noqa
        with pytest.raises(TypeError):
            assert config[1] == 1  # noqa

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
        # update自己
        _config = SwanLabConfig()
        _config.update(config)
        assert _config == config
        # update，argparse.Namespace
        _config = SwanLabConfig()
        _config.update(argparse.Namespace(a=1, b=2))
        assert _config["a"] == 1
        assert _config["b"] == 2
        # update, use kwargs
        _config = SwanLabConfig()
        _config.update(a=2, b=1)
        assert _config["a"] == 2
        assert _config["b"] == 1


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


class TestSwanLabConfigWithRun:
    """
    测试SwanLabConfig与SwanLabRun的交互
    """

    def test_use_dict(self):
        """
        正常流程，输入字典
        """
        run = SwanLabRun(
            run_config={
                "a": 1,
                "b": "mnist",
                "c/d": [1, 2, 3],
                "e/f/h": {"a": 1, "b": {"c": 2}},
            }
        )
        config = run.config
        _config = get_config()
        assert config["a"] == _config["a"] == 1
        assert config["b"] == _config["b"] == "mnist"
        assert config["c/d"] == _config["c/d"] == [1, 2, 3]

    def test_use_omegaconf(self):
        """
        正常流程，输入OmegaConf
        """
        run = SwanLabRun(
            run_config=omegaconf.OmegaConf.create(
                {
                    "a": 1,
                    "b": "mnist",
                    "c/d": [1, 2, 3],
                    "e/f/h": {"a": 1, "b": {"c": 2}},
                }
            )
        )
        config = run.config
        _config = get_config()
        assert config["a"] == _config["a"] == 1
        assert config["b"] == _config["b"] == "mnist"
        assert config["c/d"] == _config["c/d"] == [1, 2, 3]

    def test_use_argparse(self):
        """
        正常流程，输入argparse.Namespace
        """
        run = SwanLabRun(run_config=argparse.Namespace(a=1, b="mnist", c=[1, 2, 3], d={"a": 1, "b": {"c": 2}}))
        config = run.config
        _config = get_config()
        assert config["a"] == _config["a"] == 1
        assert config["b"] == _config["b"] == "mnist"
        assert config["c"] == _config["c"] == [1, 2, 3]

    def test_use_config(self):
        """
        正常流程，输入SwanLabConfig
        """
        run = SwanLabRun(
            run_config=SwanLabConfig(
                {
                    "a": 1,
                    "b": "mnist",
                    "c": [1, 2, 3],
                    "e/f/h": {"a": 1, "b": {"c": 2}},
                }
            )
        )
        config = run.config
        _config = get_config()
        assert config["a"] == _config["a"] == 1
        assert config["b"] == _config["b"] == "mnist"
        assert config["c"] == _config["c"] == [1, 2, 3]

    def test_after_finish(self):
        """
        测试在finish之后config的变化
        """
        run = SwanLabRun(
            run_config={
                "a": 1,
                "b": "mnist",
                "c/d": [1, 2, 3],
                "e/f/h": {"a": 1, "b": {"c": 2}},
            }
        )
        run.finish()
        config = run.config
        _config = get_config()
        assert len(config) == 4
        assert len(_config) == 0

    def test_error_config_input(self):
        """
        测试错误的输入
        """
        with pytest.raises(TypeError):
            SwanLabRun(run_config=1)
        with pytest.raises(TypeError):
            SwanLabRun(run_config="1")
        with pytest.raises(TypeError):
            SwanLabRun(run_config=[1, 2, 3])
        with pytest.raises(TypeError):
            SwanLabRun(run_config=(1, 2, 3))
        with pytest.raises(TypeError):
            SwanLabRun(run_config=True)
