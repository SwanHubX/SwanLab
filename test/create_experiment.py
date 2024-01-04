#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-30 15:54:55
@File: test/create_experiment.py
@IDE: vscode
@Description:
    开启一个实验
"""
import random
import swanlab as sw
import time
import numpy as np
from PIL import Image as PILImage
import os


class Enlarge1000(sw.data.BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)
        return int(self.value * 1000)

    def get_config(self, *args, **kwargs) -> dict:
        return {"color": "red"}

    def get_namespace(self, *args, **kwargs) -> str:
        return "custom"

    def get_chart_type(self) -> str:
        return self.chart.image


class Image(sw.data.BaseType):
    def get_data(self):
        print("step {}, 获取data".format(self.step))
        print(self.step, self.tag, self.settings.static_dir)

        self.image = self.preprocess(self.value)
        save_path = os.path.join(self.settings.static_dir, f"image-{self.step}.png")

        self.save(save_path)
        print("save_path:", save_path)

        save_relative_path = self.extract_path_layers(save_path)
        print("save_relative_path", save_relative_path)

        return save_relative_path

    def preprocess(self, data):
        """将不同类型的输入转换为PIL图像"""
        if isinstance(data, str):
            # 如果输入为字符串
            image = self._load_image_from_path(data)
        elif isinstance(data, np.ndarray):
            # 如果输入为numpy array
            image = self._convert_numpy_array_to_image(data)
        elif isinstance(data, PILImage.Image):
            # 如果输入为PIL.Image
            image = data
        else:
            # 以上都不是，则报错
            raise TypeError("Unsupported image type. Please provide a valid path, numpy array, or PIL.Image.")

        # 对数据做通道转换
        image = self._convert_channels(image)
        image = self.resize(image)

        return image

    def load_image_from_path(self, path):
        """判断字符串是否为正确的图像路径，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.open(path)
        except Exception as e:
            raise ValueError(f"Invalid image path: {path}") from e

    def convert_numpy_array_to_image(self, array):
        """判断np array对象是否能转换为PIL.Image，如果是则返回PIL.Image类型对象，如果不是则报错"""
        try:
            return PILImage.fromarray(array)
        except Exception as e:
            raise ValueError("Invalid numpy array for the image") from e

    def convert_channels(self, image):
        """将1通道和4通道图像转换为3通道的RGB图像。"""
        if image.mode == "L":  # 1-channel
            return image.convert("RGB")
        elif image.mode == "RGBA" or image.mode == "CMYK":  # 4-channel
            return image.convert("RGB")
        return image

    def resize(self, MAX_DIMENSION=600):
        """将图像调整大小, 保证最大边长不超过MAX_DIMENSION"""
        if max(self.image.size) > MAX_DIMENSION:
            self.image.thumbnail((MAX_DIMENSION, MAX_DIMENSION))

    def save(self, save_path):
        """将图像保存到指定路径"""
        try:
            self.image.save(save_path)
            return save_path
        except Exception as e:
            raise ValueError(f"Could not save the image to the path: {save_path}") from e

    def extract_path_layers(absolute_path: str) -> str:
        """获取绝对路径的最后三个层级的部分"""
        parts = absolute_path.split("/")
        last_three_layers = "/".join(parts[-3:])

        return last_three_layers

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Image"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.image


# 迭代次数
epochs = 10
# 学习率
lr = 0.01
# 随机偏移量
offset = random.random() / 5
# 创建一个实验
sw.init(
    description="this is a test experiment",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
    log_level="debug",
)

print("start training")

print("")

# 模拟训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    # if epoch < 10:
    #     sw.log({"loss": Enlarge1000(loss), "accuracy": acc}, step=1)
    # else:
    # sw.log({"loss": sw.data.Image("./test_image.jpg"), "accuracy": acc})

    sw.log({"loss": Image(loss), "accuracy": acc})
    # sw.log({"loss": Enlarge1000(loss), "accuracy": acc})

    # sw.log({"accuracy2": f"{acc}", "test/loss2": f"is {loss}"}, step=epochs - epoch)
    # sw.log({"loss3": loss, "accuracy3": acc}, step=1)
    # sw.log({"loss4": loss, "accuracy4": acc}, step=epoch * 2)
    time.sleep(0.5)

print("")
print("")
print("finish training")
