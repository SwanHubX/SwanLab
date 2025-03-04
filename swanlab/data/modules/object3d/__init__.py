"""SwanLab Object3D Module

This module provides classes for handling 3D data types like point clouds and meshes.

Classes:
    Object3D: Main dispatcher class for handling different types of 3D data
    PointCloud: Class for handling point cloud data
    Model3D: Class for handling 3D models

Examples:
    >>> import numpy as np
    >>> from swanlab.data.modules.object3d import Object3D
    >>> points = np.random.rand(100, 3)
    >>> obj = Object3D(points)
"""

from .model3d import Model3D
from .object3d import Object3D
from .point_cloud import PointCloud

__all__ = [
    'Object3D',
    'PointCloud',
    'Model3D',
]
