"""SwanLab Point Cloud Module

This module provides classes for handling point cloud data in various formats.

Classes:
    PointCloud: Class for handling point cloud data with different formats (XYZ, XYZC, XYZRGB)

Examples:
    Create from XYZ coordinates:
    >>> import numpy as np
    >>> points = np.random.rand(100, 3)  # 100 points with XYZ coordinates
    >>> pc = PointCloud.from_xyz(points)  # Default green color
    >>> pc = PointCloud.from_xyz(points, caption="My Points")

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
    ...     caption="Loaded Points"
    ... )
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, TypedDict

from swanlab.data.namer import hex_to_rgb, light_colors
from swanlab.toolkit import DataSuite as D, MediaBuffer, MediaType

try:
    import numpy as np

    ndarray = np.ndarray
except ImportError:
    np = None
    ndarray = None


class Box(TypedDict):
    """
    Represents a 3D bounding box with associated metadata.

    Attributes:
        color (Tuple[int, int, int]): The RGB (0..255) color of the bounding box (e.g., [255, 255, 0] for yellow).
        corners (List[List[float, float, float]]): A list of 8 corner points defining the 3D bounding box.
            Each corner is a list of three floats representing the (x, y, z) coordinates.

            The order of the corners is as follows (assuming a right-handed coordinate system):

            *   0: Bottom-Left-Back
            *   1: Bottom-Right-Back
            *   2: Bottom-Right-Front
            *   3: Bottom-Left-Front
            *   4: Top-Left-Back
            *   5: Top-Right-Back
            *   6: Top-Right-Front
            *   7: Top-Left-Front

        label (str): A string representing the object category or class label (e.g., "match").
        score (float): A confidence score representing the certainty of the bounding box's object detection.
    """

    color: Tuple[int, int, int]
    corners: List[Tuple[float, float, float]]
    label: str
    score: Optional[float]


@dataclass()
class PointCloud(MediaType):
    """Point cloud data representation.

    Supports various formats including XYZ, XYZC (XYZ + category), and XYZRGB.
    Internally stores data in XYZRGB format.

    Attributes:
        points: Point cloud data array with shape (N, 6) containing XYZRGB values
        caption: Optional description text

    Examples:
        >>> import numpy as np
        >>> points = np.random.rand(100, 3)  # XYZ format
        >>> pc = PointCloud.from_xyz(points)  # Create with default green color
        >>> pc = PointCloud.from_xyz(points, caption="My Points")  # With caption
    """

    points: ndarray  # format xyzrgb
    boxes: List[Box] = field(default_factory=list)
    step: Optional[int] = None
    caption: Optional[str] = None
    key: Optional[str] = None
    _VERSION: str = "0.1"

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
    def from_xyz(cls, points: ndarray, *, caption: Optional[str] = None, **kwargs) -> "PointCloud":
        """Create PointCloud from XYZ coordinates.

        Args:
            points: numpy array with shape (N, 3) containing XYZ coordinates
            caption: Optional description text

        Returns:
            PointCloud object with default green color

        Examples:
            >>> points = np.random.rand(100, 3)
            >>> pc = PointCloud.from_xyz(points)  # Default green
            >>> pc = PointCloud.from_xyz(points, caption="Green Points")
        """
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("XYZ array must have shape (N, 3)")

        default_color = np.array([0, 255, 0])  # green

        xyzrgb = np.zeros((points.shape[0], 6))
        xyzrgb[:, :3] = points  # copy XYZ coordinates
        xyzrgb[:, 3:] = default_color  # set default RGB (0, 255, 0) to green

        return cls(xyzrgb, caption=caption, **kwargs)

    @classmethod
    def from_xyzc(cls, points: ndarray, *, caption: Optional[str] = None, **kwargs) -> "PointCloud":
        """Create PointCloud from XYZC format (XYZ coordinates + category).

        Args:
            points: numpy array with shape (N, 4) containing XYZC values
                   where C is category index (integer)
            caption: Optional description text

        Returns:
            PointCloud object with colors mapped from categories

        Examples:
            >>> points = np.zeros((100, 4))
            >>> points[:, :3] = coordinates  # XYZ coordinates
            >>> points[:, 3] = categories    # Category labels (0,1,2...)
            >>> pc = PointCloud.from_xyzc(points)
            >>> pc = PointCloud.from_xyzc(points, caption="Segmented Points")
        """
        if points.ndim != 2 or points.shape[1] != 4:
            raise ValueError("XYZC array must have shape (N, 4)")

        # For xyzc format, map categories to predefined colors
        categories = points[:, 3].astype(int)
        colors = np.array([hex_to_rgb(c) for c in light_colors])

        xyzrgb = np.zeros((points.shape[0], 6))
        xyzrgb[:, :3] = points[:, :3]  # copy XYZ coordinates
        xyzrgb[:, 3:] = colors[categories]

        return cls(xyzrgb, caption=caption, **kwargs)

    @classmethod
    def from_xyzrgb(cls, points: ndarray, *, caption: Optional[str] = None, **kwargs) -> "PointCloud":
        """Create PointCloud from XYZRGB format.

        Args:
            points: numpy array with shape (N, 6) containing XYZRGB values
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

        return cls(points, caption=caption, **kwargs)

    @classmethod
    def from_swanlab_pts(cls, data: Dict, *, caption: Optional[str] = None, **kwargs) -> "PointCloud":
        """Create PointCloud from SwanLab pts data dictionary.

        Args:
            data: A dictionary containing 'points' (required) and optionally 'boxes'.
                  'points' can be a list of lists (XYZ, XYZC, or XYZRGB) or a NumPy array.
                  'boxes' is an optional list of dictionaries, each representing a bounding box.
            caption: Optional description text

        Returns:
            PointCloud object
        """
        if not isinstance(data, dict) or "points" not in data:
            raise ValueError("Invalid data format")

        points = data["points"]
        if isinstance(points, list):
            points = np.array(points)  # Convert list to NumPy array
        elif not isinstance(points, ndarray):
            raise TypeError("data['points'] must be a list or a NumPy array")

        if points.ndim != 2:
            raise ValueError("data['points'] must be a 2D array")

        handler = {3: cls.from_xyz, 4: cls.from_xyzc, 6: cls.from_xyzrgb}

        try:
            pc = handler[points.shape[1]](points, caption=caption, **kwargs)
        except KeyError as err:
            raise ValueError("data['points'] must have shape (N, 3), (N, 4), or (N, 6)") from err

        # Add boxes
        boxes_data = data.get("boxes")
        if boxes_data is None:
            return pc

        if not isinstance(boxes_data, list):
            raise TypeError("data['boxes'] must be a list")

        for box_data in boxes_data:
            try:
                box: Box = {
                    "color": tuple(box_data["color"]),
                    "corners": [tuple(p) for p in box_data["corners"]],
                    "label": box_data["label"],
                }
                if "score" in box_data:
                    box["score"] = box_data["score"]
                pc.append_box(box)
            except (KeyError, TypeError) as err:
                raise ValueError(f"Invalid box format: {err}") from err

        return pc

    @classmethod
    def from_swanlab_pts_json_file(cls, path: Path, *, caption: Optional[str] = None, **kwargs) -> "PointCloud":
        """Create PointCloud from SwanLab pts.json file.

        Args:
            path: Path to the .swanlab.pts.json file
            caption: Optional description text

        Returns:
            PointCloud object

        Examples:
            >>> pc = PointCloud.from_swanlab_pts_json_file(
            ...     Path("points.swanlab.pts.json"),
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

            return cls.from_swanlab_pts(points_list, caption=caption, **kwargs)

        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON format in file: {path}") from e
        except Exception as e:
            raise ValueError(f"Error reading file {path}: {str(e)}") from e

    def append_box(self, box: Box):
        """Append a bounding box to the point cloud.

        Args:
            box: The bounding box to append.
        """
        self.boxes.append(box)

    def extend_boxes(self, boxes: List[Box]):
        """Extend the list of bounding boxes with a list of boxes.

        Args:
            boxes: The list of bounding boxes to extend.
        """
        self.boxes.extend(boxes)

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
        swanlab_pts = {
            "version": self._VERSION,
            "points": points_list,
            "boxes": self.boxes,
        }
        json_str = json.dumps(swanlab_pts)
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
