"""
@author: cunyue
@file: 3d_object.py
@time: 2025/3/13 13:30
@description: 测试上传3d对象

NOTE 你需要下载下面的文件放到当前文件目录的assets文件夹下，才能运行这个测试
- 3D点云数据（带标注）: https://drive.google.com/file/d/1mFill-BXw3cirPHwIHndb1wNX4pWvSXb/view
"""

import json
import os

# noinspection PyPackageRequirements
import numpy as np

import swanlab

current_path = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------- 点云数据/swanlab点云文件路径 ----------------------------------

point_cloud_path = os.path.join(current_path, "assets/data.swanlab.pts.json")

with open(point_cloud_path, "r") as f:
    point_cloud = json.load(f)


# ---------------------------------- 从代码生成的点云矩阵 ----------------------------------


def mock_object_3d_array() -> np.ndarray:
    return np.random.rand(100, 3)


# ---------------------------------- 实验 ----------------------------------


swanlab.init(project="3d-object-project", config={"epochs": 5}, mode="cloud")


for epoch in range(2, swanlab.config.epochs):
    swanlab.log(
        {
            "code/random": [swanlab.Object3D(mock_object_3d_array()) for _ in range(epoch)],
            "epoch": epoch,
            "code/point_cloud": swanlab.Object3D(point_cloud),
        }
    )

swanlab.finish()
