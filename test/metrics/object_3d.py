"""
@author: cunyue
@file: 3d_object.py
@time: 2025/3/13 13:30
@description: 测试上传3d对象
"""

# noinspection PyPackageRequirements
import numpy as np

import swanlab


def mock_object_3d_array() -> np.ndarray:
    return np.random.rand(100, 3)


swanlab.init(project="3d-object-project", config={"epochs": 5}, mode="cloud")


for epoch in range(2, swanlab.config.epochs):
    swanlab.log({"object-3d": [swanlab.Object3D(mock_object_3d_array()) for _ in range(epoch)]})

swanlab.finish()
