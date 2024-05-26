#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/26 16:03
@File: pytest_main.py
@IDE: pycharm
@Description:
    测试SwanLabRun主类
"""
from swanlab.data.run.main import SwanLabRun, get_run, SwanLabRunState, swanlog
from swanlab import Image, Audio, Text
from nanoid import generate
from tutils import clear, TEMP_PATH
from PIL import Image as PILImage
import torch
import soundfile as sf
import numpy as np
import pytest
import random
import os


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


class TestSwanLabRunInit:

    def test_before_init(self):
        run = get_run()
        assert run is None
        assert SwanLabRun.get_state() == SwanLabRunState.NOT_STARTED

    def test_after_init(self):
        run = SwanLabRun(generate())
        assert swanlog.installed is False
        assert run is not None
        assert get_run().__str__() == run.__str__()
        _run = run.finish()
        assert get_run() is None
        assert _run.__str__() == run.__str__()
        assert SwanLabRunState.SUCCESS == _run.state

    def test_duplicate_init(self):
        run = SwanLabRun(generate())
        with pytest.raises(RuntimeError) as e:
            SwanLabRun(generate())
        assert swanlog.installed is False
        assert str(e.value) == "SwanLabRun has been initialized"
        assert run.__str__() == get_run().__str__()
        _run = run.finish()
        assert swanlog.installed is False


class TestSwanLabRunLog:
    """
    测试SwanLabRun的日志解析功能
    1. 输入为字典
    2. 可包含参数step
    3. 输入的字典的key必须为字符串
    4. 输出解析后的数据对象，包含额外的信息，step等信息
    5. 如果某一个解析失败，对应的key存在，但返回的数据为None
    """

    @pytest.mark.parametrize(
        "data",
        [
            random.randint(1, 100),
            -random.randint(1, 100),
            random.random(),
        ],
    )
    def test_log_number_ok(self, data):
        """
        测试解析一个正常的数字
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": data})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert len(a.metric) == 3
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        assert a.metric["data"] == data
        assert a.step == 0
        assert a.epoch == 1
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_nan(self):
        """
        测试解析一个nan
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": float("nan")})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert a.metric is None
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error["data_class"] == "NaN"
        assert ac.error["excepted"] == ["float", "int"]
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_str(self):
        """
        测试解析其他字符串
        """
        """
        测试解析一个nan
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": "abc"})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert a.metric is None
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error["data_class"] == "str"
        assert ac.error["excepted"] == ["float", "int"]
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_use_step(self):
        """
        测试解析一个数字，使用step
        """
        run = SwanLabRun()
        metric_dict = run.log({"a": 1}, step=1)
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert len(a.metric) == 3
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        assert a.metric["data"] == 1
        assert a.step == 1
        assert a.epoch == 1
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "default"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "default"
        assert ac.namespace == "default"
        assert ac.reference == "step"

    def test_log_number_use_step_duplicate(self):
        """
        测试解析一个数字，使用step，但是重复执行
        """
        run = SwanLabRun()
        run.log({"a": 1}, step=1)
        metric_dict = run.log({"a": 1}, step=1)
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        assert a.error is True

    def test_log_number_use_prefix(self):
        """
        测试解析一个数字，使用prefix
        """
        run = SwanLabRun()
        prefix_1 = generate()
        prefix_2 = generate() + "/" + generate()
        key1 = f"{prefix_1}/a"
        key2 = f"{prefix_2}/a"
        metric_dict = run.log({key1: 1, key2: 1})
        assert len(metric_dict) == 2
        for key in metric_dict:
            assert metric_dict[key] is not None
            a = metric_dict[key]
            ac = a.column_info
            # ---------------------------------- 指标信息 ----------------------------------
            assert len(a.metric) == 3
            assert "index" in a.metric
            assert "create_time" in a.metric
            assert "data" in a.metric
            assert a.metric["data"] == 1
            assert a.step == 0
            assert a.epoch == 1
            # ---------------------------------- 列信息 ----------------------------------
            assert ac.data_type == "default"
            assert ac.error is None
            # 默认排在最前面
            assert ac.sort is None
            assert ac.config == {}
            assert ac.key == key
            assert ac.chart_type == "default"
            assert ac.namespace == key.split("/")[0]
            assert ac.reference == "step"

    def test_log_image_path(self):
        """
        测试解析一个图片, 使用文件路径
        """
        # 创建随机图像，并保存到TEMP目录
        random_im = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
        random_im_pil = PILImage.fromarray(random_im)
        save_path = os.path.join(TEMP_PATH, "a.jpg")
        random_im_pil.save(save_path)

        run = SwanLabRun()
        metric_dict = run.log({"a": Image(save_path)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_image_np(self):
        """
        测试解析一个图片，使用np.ndarray
        """

        random_im = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)

        run = SwanLabRun()
        metric_dict = run.log({"a": Image(random_im)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_image_pil(self):
        """
        测试解析一个图片，使用PIL.Image
        """

        random_im = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)
        random_im_pil = PILImage.fromarray(random_im)

        run = SwanLabRun()
        metric_dict = run.log({"a": Image(random_im_pil)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_image_plt(self):
        """
        测试解析一个图片，使用Matplotlib
        """
        import matplotlib.pyplot as plt

        x = [1, 2, 3, 4, 5]
        y = [2, 3, 5, 7, 11]
        plt.plot(x, y)
        plt.title("Examples")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")

        run = SwanLabRun()
        metric_dict = run.log({"a": Image(plt)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_image_tensor(self):
        """
        测试解析一个图片，使用pytorch tensor
        """
        random_im = torch.randn(4, 3, 64, 64)

        run = SwanLabRun()
        metric_dict = run.log({"a": Image(random_im)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_image_batch(self):
        """
        测试解析一次log一批图片
        """
        random_im = np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8)

        run = SwanLabRun()
        metric_dict = run.log({"a": [Image(random_im)] * 3})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "image"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "image"
        assert ac.namespace == "Image"
        assert ac.reference == "step"

    def test_log_audio_path(self):
        """
        测试解析一个音频，使用文件路径
        """
        sample_rate = 44100
        t = np.linspace(0, 1, sample_rate, False)
        frequency = 440
        audio_signal = 0.5 * np.sin(2 * np.pi * frequency * t)
        save_path = os.path.join(TEMP_PATH, "output.wav")
        sf.write(save_path, audio_signal, sample_rate)

        run = SwanLabRun()
        metric_dict = run.log({"a": Audio(save_path)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "audio"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "audio"
        assert ac.namespace == "Audio"
        assert ac.reference == "step"

    def test_log_audio_np(self):
        """
        测试解析一个音频，使用np.ndarray
        """
        random_audio = np.random.randn(2, 100000)

        run = SwanLabRun()
        metric_dict = run.log({"a": Audio(random_audio)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "audio"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "audio"
        assert ac.namespace == "Audio"
        assert ac.reference == "step"

    def test_log_audio_batch(self):
        """
        测试解析一批音频
        """
        random_audio = np.random.randn(2, 100000)

        run = SwanLabRun()
        metric_dict = run.log({"a": [Audio(random_audio)] * 3})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "audio"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "audio"
        assert ac.namespace == "Audio"
        assert ac.reference == "step"

    def test_log_text_str(self):
        """
        测试解析单个字符串文本
        """
        random_text = "Hello World"

        run = SwanLabRun()
        metric_dict = run.log({"a": Text(random_text)})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "text"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "text"
        assert ac.namespace == "Text"
        assert ac.reference == "step"

    def test_log_text_bacth(self):
        """
        测试解析一批文本
        """
        random_text = "Hello World"

        run = SwanLabRun()
        metric_dict = run.log({"a": [Text(random_text)] * 3})
        assert metric_dict["a"] is not None
        assert len(metric_dict) == 1
        a = metric_dict["a"]
        ac = a.column_info
        # ---------------------------------- 指标信息 ----------------------------------
        assert "index" in a.metric
        assert "create_time" in a.metric
        assert "data" in a.metric
        # ---------------------------------- 列信息 ----------------------------------
        assert ac.data_type == "text"
        assert ac.error is None
        # 默认排在最前面
        assert ac.sort == 0
        assert ac.config == {}
        assert ac.key == "a"
        assert ac.chart_type == "text"
        assert ac.namespace == "Text"
        assert ac.reference == "step"


class TestSwanLabRunState:
    """
    测试SwanLabRun的状态变化
    """

    def test_not_started(self):
        assert SwanLabRun.get_state() == SwanLabRunState.NOT_STARTED

    def test_running(self):
        run = SwanLabRun()
        assert run.state == SwanLabRunState.RUNNING
        assert run.is_running is True

    def test_crashed(self):
        run = SwanLabRun()
        run.finish(SwanLabRunState.CRASHED, error="error")
        assert run.state == SwanLabRunState.CRASHED
        assert run.is_crashed is True

    def test_success(self):
        run = SwanLabRun()
        run.finish(SwanLabRunState.SUCCESS)
        assert run.state == SwanLabRunState.SUCCESS
        assert run.is_success is True
