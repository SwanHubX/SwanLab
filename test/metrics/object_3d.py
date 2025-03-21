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

cloud_point_path = os.path.join(current_path, "assets/data.swanlab.pts.json")

with open(cloud_point_path, "r") as f:
    cloud_point = json.load(f)


# ---------------------------------- 从代码生成的点云矩阵 ----------------------------------


def mock_object_3d_array() -> np.ndarray:
    return np.random.rand(100, 3)


# ---------------------------------- 实验 ----------------------------------


swanlab.init(project="3d-object-project", config={"epochs": 10}, mode="cloud")


for epoch in range(2, swanlab.config.epochs):
    swanlab.log(
        {
            "code/random": [swanlab.Object3D(mock_object_3d_array()) for _ in range(epoch)],
            "epoch": epoch,
            "code/point_cloud": swanlab.Object3D.from_point_data(
                points=cloud_point['points'], boxes=cloud_point['boxes'], caption="Point Cloud from cloud"
            ),
            "file/point_cloud": swanlab.Object3D(cloud_point_path, caption="Point Cloud from file"),
        }
    )
