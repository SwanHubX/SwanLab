from swankit.core import MediaType, MediaBuffer, DataSuite as D
from typing import Any, override, Optional, Dict, Tuple
from functools import cached_property
from swanlab.data.run.namer import hex_to_rgb, light_colors
import json

try:
    import numpy as np
    PointsArray = np.ndarray
except ImportError:
    np = None
    PointsArray = Any


class PointCloud(MediaType):
    """A class representing 3D point cloud data.

    Supports three formats:
    - xyz: coordinates only
    - xyzc: coordinates with category
    - xyzrgb: coordinates with RGB colors

    Args:
        points: numpy array of point cloud data
        caption: optional description of the point cloud
    """

    def __init__(self, points: PointsArray,  caption: Optional[str] = None):
        # 在使用时检查numpy是否已安装
        if np is None:
            raise ImportError(
                "Numpy is required for PointCloud class. "
                "Please install it with: pip install numpy."
            )

        super().__init__()

        self.points = PointsData(points).convert_to_xyzrgb()
        self.caption = D.check_caption(caption)

        self.buffer = MediaBuffer()
        points_list = self.points.take().tolist()
        json.dump(points_list, self.buffer)

    @classmethod
    def from_file(cls, file_path: str, caption: Optional[str] = None) -> "PointCloud":
        """从JSON文件创建点云对象

        Args:
            file_path: JSON文件路径
            caption: 可选的点云描述

        Returns:
            PointCloud 对象

        Raises:
            ImportError: 如果numpy未安装
            ValueError: 如果文件格式不正确或数据无效
            FileNotFoundError: 如果文件不存在
            """
        if np is None:
            raise ImportError(
                "Numpy is required for PointCloud class. "
                "Please install it with: pip install numpy."
            )

        try:
            with open(file_path, 'r') as f:
                points_list = json.load(f)

            # 检查是否为嵌套列表
            if not isinstance(points_list, list) or not all(isinstance(p, list) for p in points_list):
                raise ValueError(
                    "JSON file must contain a list of point lists")

            # 转换为numpy数组
            points = np.array(points_list)

        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find file: {file_path}")
        except json.JSONDecodeError:
            raise ValueError(f"Invalid JSON format in file: {file_path}")
        except Exception as e:
            raise ValueError(f"Error loading point cloud data: {str(e)}")

        # 通过创建实例时的验证来检查数据格式
        return cls(points, caption)

    @override
    def parse(self) -> Tuple[str, MediaBuffer]:
        """生成文件名并返回数据"""
        # 根据数据生成哈希值作为文件名的一部分
        hash_name = D.get_hash_by_ndarray(self.points.take())[:16]
        save_name = f"pointscloud-step{self.step}-{hash_name}.json"
        return save_name, self.buffer  # 返回文件名和数据

    @override
    def get_more(self) -> Optional[Dict[str, str]]:
        """返回额外信息"""
        return {"caption": self.caption} if self.caption else None

    @override
    def get_section(self) -> str:
        """返回分组名"""
        return "PointsCloud"

    @override
    def get_chart(self):
        """返回图表类型"""
        return self.Chart.ChartItem


class PointsData:
    POINT_AVAILABLE_FORMAT = {3: 'xyz', 4: 'xyzc', 6: 'xyzrgb'}
    DEFAULT_COLOR = (255, 255, 255)

    def __init__(self, points: PointsArray):
        self.validate_points(points)
        self.points = points

    def take(self) -> PointsArray:
        return self.points

    def convert_to_xyzrgb(self) -> "PointsData":
        """将点云数据转换为xyzrgb格式
        对于xyz格式: 添加默认白色(255, 255, 255)
        对于xyzc格式: 根据类别映射到预定义的颜色
        对于xyzrgb格式: 保持不变
        """
        if self.format == 'xyzrgb':
            return self

        # 创建新的数组
        n_points = len(self.points)
        xyzrgb = np.zeros((n_points, 6), dtype=self.points.dtype)
        # 复制xyz坐标
        xyzrgb[:, :3] = self.points[:, :3]

        if self.format == 'xyz':
            # 使用默认灰色
            xyzrgb[:, 3:] = PointsData.DEFAULT_COLOR
            return PointsData(xyzrgb)

        # xyzc格式 预定义16种颜色 根据类别获取颜色
        categories = self.points[:, 3].astype(int)
        colors = np.array([hex_to_rgb(c) for c in light_colors])
        xyzrgb[:, 3:] = colors[categories]
        return PointsData(xyzrgb)

    @cached_property
    def format(self) -> str:
        return PointsData.POINT_AVAILABLE_FORMAT[self.points.shape[1]]

    @classmethod
    def from_xyz(cls, points: PointsArray) -> "PointsData":
        """从xyz数据构建"""
        p = cls(points)
        if p.format != "xyz":
            raise ValueError("XYZ data must have 3 features")
        return p

    @classmethod
    def from_xyzc(cls, points: PointsArray) -> "PointsData":
        """从xyzc数据构建"""
        p = cls(points)
        if p.format != "xyzc":
            raise ValueError("XYZC data must have 4 features")
        return p

    @classmethod
    def from_xyzrgb(cls, points: PointsArray) -> "PointsData":
        """从xyzrgb数据构建"""
        p = cls(points)
        if p.format != "xyzrgb":
            raise ValueError("XYZRGB data must have 6 features")
        return p

    @staticmethod
    def validate_points(points: PointsArray) -> None:
        if points is None:
            raise ValueError("Points array cannot be None")

        if len(points) == 0:
            raise ValueError("Points array cannot be empty")

        # 检查维度
        if points.ndim != 2:
            raise ValueError(
                f"Points must be a 2D array, got {points.ndim}D array."
            )

        # 检查输入点是否是所需要的格式
        n_features = points.shape[1]
        if n_features not in PointsData.POINT_AVAILABLE_FORMAT.keys():
            raise ValueError(
                "Point cloud data must be in one of these formats:\n"
                "- Nx3: [x, y, z] coordinates only\n"
                "- Nx4: [x, y, z, category] with category in range int[0,15]\n"
                "- Nx6: [x, y, z, r, g, b] with RGB values in range int[0,255]\n"
                f"But got {n_features} shaped point."
            )

        # 检查是否是数字
        if not np.issubdtype(points.dtype, np.number):
            raise ValueError("points must contain numeric values")

        # 检查 category 是否在 int[0, 15]
        if n_features == 4:
            categories = points[:, 3]
            # 检查是否为整数
            if not np.all(categories == categories.astype(int)):
                raise ValueError("Categories must be integers")

            # 检查范围
            if not np.all((categories >= 0) & (categories <= 15)):
                raise ValueError("Categories must be in range [0, 15]")

        # 检查 rgb 是否在 int[0, 255]
        elif n_features == 6:
            rgb = points[:, 3:6]
            # 检查RGB是否合法
            if not np.all(rgb == rgb.astype(int)):
                raise ValueError("RGB values must be integers")

            # 检查范围
            if not np.all((rgb >= 0) & (rgb <= 255)):
                raise ValueError("RGB values must be in range [0, 255]")
