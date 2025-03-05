"""SwanLab Point Cloud Module

This module provides classes for handling point cloud data in various formats.

Classes:
    PointCloud: Class for handling point cloud data with different formats (XYZ, XYZC, XYZRGB)

Examples:
    Create from XYZ coordinates:
    >>> import numpy as np
    >>> points = np.random.rand(100, 3)  # 100 points with XYZ coordinates
    >>> pc = PointCloud.from_xyz(points)  # Default green color
    >>> pc = PointCloud.from_xyz(points, step=1, caption="My Points")

    Create from XYZC format (XYZ + category):
    >>> points_with_category = np.zeros((100, 4))  # XYZ + category
    >>> points_with_category[:, :3] = points  # XYZ coordinates
    >>> points_with_category[:, 3] = categories  # Category labels (0,1,2...)
    >>> pc = PointCloud.from_xyzc(points_with_category)  # Auto color mapping

    Create from XYZRGB format:
    >>> points_with_color = np.zeros((100, 6))  # XYZRGB format
    >>> points_with_color[:, :3] = points  # XYZ coordinates
    >>> points_with_color[:, 3:] = colors  # RGB values (0-255)
    >>> pc = PointCloud.from_xyzrgb(points_with_color)

    Load from SwanLab pts.json file:
    >>> from pathlib import Path
    >>> pc = PointCloud.from_swanlab_pts_json_file(
    ...     Path("./file/points.swanlab.pts.json"),
    ...     step=1,
    ...     caption="Loaded Points"
    ... )
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

from swankit.core.data import DataSuite as D
from swankit.core.data import MediaBuffer, MediaType

from swanlab.data.namer import hex_to_rgb, light_colors

try:
    import numpy as np
except ImportError:
    np = None


@dataclass(frozen=True)
class PointCloud(MediaType):
    """Point cloud data representation.

    Supports various formats including XYZ, XYZC (XYZ + category), and XYZRGB.
    Internally stores data in XYZRGB format.

    Attributes:
        points: Point cloud data array with shape (N, 6) containing XYZRGB values
        step: Optional step number for visualization sequences
        caption: Optional description text

    Examples:
        >>> import numpy as np
        >>> points = np.random.rand(100, 3)  # XYZ format
        >>> pc = PointCloud.from_xyz(points)  # Create with default green color
        >>> pc = PointCloud.from_xyz(points, step=1)  # With step number
        >>> pc = PointCloud.from_xyz(points, caption="My Points")  # With caption
    """

    points: np.ndarray  # format xyzrgb
    step: Optional[int] = None
    caption: Optional[str] = None

    def __post_init__(self):
        """Validate input data after initialization"""
        self._check_numpy()

        if self.points.ndim != 2 or self.points.shape[1] != 6:
            raise ValueError("PointCloud(points) points must have shape (N, 6)")

    @staticmethod
    def _check_numpy() -> None:
        if np is None:
            raise ImportError("Numpy is required for PointCloud class. Please install it with: pip install numpy")

    @classmethod
    def from_xyz(
        cls, points: np.ndarray, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "PointCloud":
        """Create PointCloud from XYZ coordinates.

        Args:
            points: numpy array with shape (N, 3) containing XYZ coordinates
            step: Optional step number for visualization sequences
            caption: Optional description text

        Returns:
            PointCloud object with default green color

        Examples:
            >>> points = np.random.rand(100, 3)
            >>> pc = PointCloud.from_xyz(points)  # Default green
            >>> pc = PointCloud.from_xyz(points, step=1, caption="Green Points")
        """
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("XYZ array must have shape (N, 3)")

        default_color = np.array([0, 255, 0])  # green

        xyzrgb = np.zeros((points.shape[0], 6))
        xyzrgb[:, :3] = points  # copy XYZ coordinates
        xyzrgb[:, 3:] = default_color  # set default RGB (0, 255, 0) to green

        return cls(xyzrgb, step=step, caption=caption, **kwargs)

    @classmethod
    def from_xyzc(
        cls, points: np.ndarray, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "PointCloud":
        """Create PointCloud from XYZC format (XYZ coordinates + category).

        Args:
            points: numpy array with shape (N, 4) containing XYZC values
                   where C is category index (integer)
            step: Optional step number for visualization sequences
            caption: Optional description text

        Returns:
            PointCloud object with colors mapped from categories

        Examples:
            >>> points = np.zeros((100, 4))
            >>> points[:, :3] = coordinates  # XYZ coordinates
            >>> points[:, 3] = categories    # Category labels (0,1,2...)
            >>> pc = PointCloud.from_xyzc(points)
            >>> pc = PointCloud.from_xyzc(points, step=1, caption="Segmented Points")
        """
        if points.ndim != 2 or points.shape[1] != 4:
            raise ValueError("XYZC array must have shape (N, 4)")

        # For xyzc format, map categories to predefined colors
        categories = points[:, 3].astype(int)
        colors = np.array([hex_to_rgb(c) for c in light_colors])

        xyzrgb = np.zeros((points.shape[0], 6))
        xyzrgb[:, :3] = points[:, :3]  # copy XYZ coordinates
        xyzrgb[:, 3:] = colors[categories]

        return cls(xyzrgb, step=step, caption=caption, **kwargs)

    @classmethod
    def from_xyzrgb(
        cls, points: np.ndarray, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "PointCloud":
        """Create PointCloud from XYZRGB format.

        Args:
            points: numpy array with shape (N, 6) containing XYZRGB values
            step: Optional step number for visualization sequences
            caption: Optional description text

        Returns:
            PointCloud object

        Examples:
            >>> points = np.zeros((100, 6))
            >>> points[:, :3] = coordinates  # XYZ coordinates
            >>> points[:, 3:] = colors       # RGB values (0-255)
            >>> pc = PointCloud.from_xyzrgb(points)
            >>> pc = PointCloud.from_xyzrgb(points, step=1, caption="Colored Points")
        """
        if points.ndim != 2 or points.shape[1] != 6:
            raise ValueError("XYZRGB array must have shape (N, 6)")

        return cls(points, step=step, caption=caption, **kwargs)

    @classmethod
    def from_swanlab_pts_json_file(
        cls, path: Path, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "PointCloud":
        """Create PointCloud from SwanLab pts.json file.

        The file should contain a list of points in XYZRGB format:
        [[x1,y1,z1,r1,g1,b1], [x2,y2,z2,r2,g2,b2], ...]

        Args:
            path: Path to the .swanlab.pts.json file
            step: Optional step number for visualization sequences
            caption: Optional description text

        Returns:
            PointCloud object

        Examples:
            >>> pc = PointCloud.from_swanlab_pts_json_file(
            ...     Path("points.swanlab.pts.json"),
            ...     step=1,
            ...     caption="Loaded Points"
            ... )
        """
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")

        if not path.is_file():
            raise ValueError(f"Path is not a file: {path}")

        try:
            with open(path) as f:
                points_list = json.load(f)

            if not isinstance(points_list, list) or not points_list:
                raise ValueError("Invalid file format: expected non-empty list")

            if not isinstance(points_list[0], list) or len(points_list[0]) != 6:
                raise ValueError("Invalid point format: expected [x,y,z,r,g,b]")

            points = np.array(points_list)
            return cls(points, step=step, caption=caption, **kwargs)

        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON format in file: {path}") from e
        except Exception as e:
            raise ValueError(f"Error reading file {path}: {str(e)}") from e

    # ---------------------------------- override ----------------------------------

    def parse(self) -> Tuple[str, MediaBuffer]:
        """Convert point cloud to buffer for transmission.

        Returns:
            Tuple containing:
            - File name with format: pointcloud-step{step}-{hash}.swanlab.pts.json
            - MediaBuffer containing the point cloud data
        """
        buffer = MediaBuffer()
        points_list = self.points.tolist()
        json_str = json.dumps(points_list)
        buffer.write(json_str.encode())

        hash_name = D.get_hash_by_ndarray(self.points)[:16]
        save_name = f"pointcloud-step{self.step}-{hash_name}.swanlab.pts.json"

        return save_name, buffer

    def get_chart(self) -> MediaType.Chart:
        """Return chart type for visualization"""
        return MediaType.Chart.OBJECT3D

    def get_section(self) -> str:
        """Return section name for organization"""
        return "PointsCloud"

    def get_more(self) -> Optional[Dict[str, str]]:
        """Return additional information (caption)"""
        return {"caption": self.caption} if self.caption else None
