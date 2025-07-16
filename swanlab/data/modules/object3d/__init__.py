"""SwanLab Object3D Module

This module provides classes for handling 3D data types like point clouds and meshes.

Classes:
    Object3D: Main dispatcher class for handling different types of 3D data
    PointCloud: Class for handling point cloud data with XYZ, XYZC, and XYZRGB formats
    Model3D: Class for handling 3D model files like GLB
    Molecule: Class for handling molecule data from various formats by RDKit

Examples:
    # Create point cloud from XYZ coordinates
    >>> import numpy as np
    >>> from swanlab.data.modules.object3d import Object3D
    >>> points = np.random.rand(100, 3)  # XYZ format
    >>> obj = Object3D(points)  # Returns PointCloud object

    # Create point cloud from XYZC data (coordinates + category)
    >>> points_c = np.random.rand(100, 4)  # XYZC format
    >>> obj = Object3D(points_c)  # Category colors mapped automatically

    # Create point cloud from XYZRGB data
    >>> points_rgb = np.random.rand(100, 6)  # XYZRGB format
    >>> obj = Object3D(points_rgb, step=1, caption="My Points")

    # Load from files
    >>> obj = Object3D("model.glb")  # Load 3D model
    >>> obj = Object3D("cloud.swanlab.pts.json")  # Load point cloud

    # Direct class usage
    >>> from swanlab.data.modules.object3d import PointCloud, Model3D
    >>> pc = PointCloud.from_xyz(points, step=1)
    >>> model = Model3D.from_glb_file("model.glb", caption="My Model")
"""

from .model3d import Model3D
from .molecule import Molecule
from .object3d import Object3D
from .point_cloud import PointCloud

__all__ = [
    'Object3D',
    'PointCloud',
    'Model3D',
    'Molecule',
]
