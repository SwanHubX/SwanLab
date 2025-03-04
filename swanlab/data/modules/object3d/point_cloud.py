import json
from dataclasses import dataclass
from functools import once
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
            raise ImportError("Numpy is required for PointCloud class. " "Please install it with: pip install numpy.")

    @classmethod
    def from_xyz(
        cls, points: np.ndarray, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "PointCloud":
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
        if points.ndim != 2 or points.shape[1] != 6:
            raise ValueError("XYZRGB array must have shape (N, 6)")

        return cls(points, step=step, caption=caption, **kwargs)

    # ---------------------------------- override ----------------------------------

    @once
    def parse(self) -> Tuple[str, MediaBuffer]:
        buffer = MediaBuffer()
        points_list = self.points.tolist()
        json_str = json.dumps(points_list)
        buffer.write(json_str.encode())

        hash_name = D.get_hash_by_ndarray(self.points)[:16]
        save_name = f"pointcloud-step{self.step}-{hash_name}.swanlab.pts.json"

        return save_name, buffer

    def get_chart(self) -> MediaType.Chart:
        return MediaType.Chart.OBJECT3D

    def get_section(self) -> str:
        return "PointsCloud"

    def get_more(self) -> Optional[Dict[str, str]]:
        return {"caption": self.caption} if self.caption else None
