"""
@author: cunyue
@file: __init__.py
@time: 2026/4/30
@description: 3D 对象处理模块，支持点云和3D模型
"""

import hashlib
import json
from io import BytesIO
from pathlib import Path
from typing import Optional

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.object3d import Object3DDataType

# ---------- 颜色映射常量（用于 xyzc -> xyzrgb 转换） ----------

_LIGHT_COLORS = [
    "#528d59",
    "#587ad2",
    "#c24d46",
    "#9cbe5d",
    "#6ebad3",
    "#dfb142",
    "#6d4ba4",
    "#8cc5b7",
    "#892d58",
    "#40877c",
    "#d0703c",
    "#d47694",
    "#e3b292",
    "#b15fbb",
    "#905f4a",
    "#989fa3",
]


def _hex_to_rgb(hex_color: str):
    hex_color = hex_color.strip().lstrip("#").strip()
    return int(hex_color[:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)


_LIGHT_COLORS_RGB = [_hex_to_rgb(c) for c in _LIGHT_COLORS]

_PTS_VERSION = "0.1"

# ---------- 文件扩展名 -> 类型映射 ----------

_FILE_EXT_MAP = {
    ".swanlab.pts.json": "pointcloud",
    ".glb": "model3d",
}


class Object3D(TransformMedia):
    """3D 对象转换类，支持点云和3D模型。

    支持的输入类型:
    - numpy.ndarray: 点云数据，形状 (N, 3) XYZ / (N, 4) XYZC / (N, 6) XYZRGB
    - dict: {"points": ndarray, "boxes": [...]} 格式的点云数据
    - str/Path: 文件路径，支持 .glb 和 .swanlab.pts.json
    - Object3D: 套娃加载
    """

    def __init__(self, data: Object3DDataType, caption: CaptionType = None):
        super().__init__()

        # 套娃加载
        attrs = self._unwrap(data)
        if attrs:
            self.buffer: BytesIO = attrs["buffer"]
            self.file_type: str = attrs["file_type"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        self.caption = caption

        # ---------- 按输入类型分发处理 ----------
        np = vendor.np

        if isinstance(data, dict):
            self.buffer, self.file_type = self._handle_dict(data, np)
        elif isinstance(data, (str, Path)):
            self.buffer, self.file_type = self._handle_file(Path(data))
        elif isinstance(data, np.ndarray):
            self.buffer, self.file_type = self._handle_ndarray(data, np)
        else:
            raise TypeError(f"Unsupported input type: {type(data)}. Expected numpy.ndarray, dict, str, or Path.")

    # ---------- ndarray 处理 ----------

    @staticmethod
    def _handle_ndarray(data, np) -> tuple:
        """将 numpy 点云数据转换为 swanlab.pts.json 格式的 buffer。"""
        if data.ndim != 2 or data.shape[1] not in (3, 4, 6):
            raise ValueError(
                f"Unsupported array shape: {data.shape}. "
                f"Expected 2D array with 3 (XYZ), 4 (XYZC), or 6 (XYZRGB) channels."
            )

        channels = data.shape[1]
        if channels == 3:
            xyzrgb = np.zeros((data.shape[0], 6))
            xyzrgb[:, :3] = data
            xyzrgb[:, 3:] = [0, 255, 0]  # 默认绿色
        elif channels == 4:
            categories = data[:, 3].astype(int)
            colors = np.array(_LIGHT_COLORS_RGB)
            xyzrgb = np.zeros((data.shape[0], 6))
            xyzrgb[:, :3] = data[:, :3]
            xyzrgb[:, 3:] = colors[categories % len(colors)]
        else:  # channels == 6
            xyzrgb = data

        pts_data = {
            "version": _PTS_VERSION,
            "points": xyzrgb.tolist(),
            "boxes": [],
        }
        return _json_to_buffer(pts_data), "swanlab.pts.json"

    # ---------- dict 处理 ----------

    @staticmethod
    def _handle_dict(data: dict, np) -> tuple:
        """将 {"points": ..., "boxes": [...]} 字典转换为 buffer。"""
        if not isinstance(data, dict) or "points" not in data:
            raise ValueError("Dict input must contain 'points' key.")

        points = data["points"]
        if isinstance(points, list):
            points = np.array(points)
        elif not isinstance(points, np.ndarray):
            raise TypeError("data['points'] must be a list or numpy.ndarray")

        if points.ndim != 2 or points.shape[1] not in (3, 4, 6):
            raise ValueError("data['points'] must have shape (N, 3), (N, 4), or (N, 6)")

        channels = points.shape[1]
        if channels == 3:
            xyzrgb = np.zeros((points.shape[0], 6))
            xyzrgb[:, :3] = points
            xyzrgb[:, 3:] = [0, 255, 0]
        elif channels == 4:
            categories = points[:, 3].astype(int)
            colors = np.array(_LIGHT_COLORS_RGB)
            xyzrgb = np.zeros((points.shape[0], 6))
            xyzrgb[:, :3] = points[:, :3]
            xyzrgb[:, 3:] = colors[categories % len(colors)]
        else:
            xyzrgb = points

        # 处理 boxes
        boxes_data = data.get("boxes", [])
        boxes = []
        for box in boxes_data:
            b = {
                "color": list(box["color"]),
                "corners": [list(p) for p in box["corners"]],
                "label": box["label"],
            }
            if "score" in box:
                b["score"] = box["score"]
            boxes.append(b)

        pts_data = {
            "version": _PTS_VERSION,
            "points": xyzrgb.tolist(),
            "boxes": boxes,
        }
        return _json_to_buffer(pts_data), "swanlab.pts.json"

    # ---------- 文件处理 ----------

    @staticmethod
    def _handle_file(path: Path) -> tuple:
        """根据文件扩展名加载文件到 buffer。"""
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        if not path.is_file():
            raise ValueError(f"Path is not a file: {path}")

        suffixes = path.suffixes
        if not suffixes:
            raise ValueError(f"File has no extension: {path}")

        # 渐进式匹配后缀（如 .swanlab.pts.json）
        all_tried = ["".join(suffixes[i:]).lower() for i in range(len(suffixes))]

        for suffix in all_tried:
            file_type = _FILE_EXT_MAP.get(suffix)
            if file_type is not None:
                if file_type == "pointcloud":
                    return _load_pts_json_file(path), "swanlab.pts.json"
                elif file_type == "model3d":
                    return _load_glb_file(path), "glb"

        raise ValueError(
            f"Unsupported file type: {path.name}. Supported extensions: {', '.join(sorted(_FILE_EXT_MAP.keys()))}"
        )

    # ---------- 工厂方法 ----------

    @classmethod
    def from_point_data(cls, points, boxes=None, caption: CaptionType = None) -> "Object3D":
        """从点云数据和检测框创建 Object3D 对象。

        Args:
            points: 点云数组，支持 (N,3) XYZ / (N,4) XYZC / (N,6) XYZRGB
            boxes: 可选的3D检测框列表，每个框包含 color、corners、label、score
            caption: 可选的说明文字
        """
        return cls({"points": points, "boxes": boxes or []}, caption=caption)

    # ---------- 抽象方法实现 ----------

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_OBJECT3D

    def transform(self, *, step: int, path: Path) -> MediaItem:
        content = self.buffer.getvalue()
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.{self.file_type}"
        fs.safe_write(path / filename, content, mode="wb")
        return MediaItem(filename=filename, sha256=sha256, size=len(content), caption=self.caption or "")


# ---------- 辅助函数 ----------


def _json_to_buffer(data: dict) -> BytesIO:
    """将字典序列化为 JSON 并写入 BytesIO。"""
    buf = BytesIO()
    buf.write(json.dumps(data).encode())
    buf.seek(0)
    return buf


def _load_pts_json_file(path: Path) -> BytesIO:
    """加载 .swanlab.pts.json 文件到 buffer。"""
    with open(path, "rb") as f:
        content = f.read()
    buf = BytesIO(content)
    buf.seek(0)
    return buf


def _load_glb_file(path: Path) -> BytesIO:
    """加载 .glb 文件到 buffer。"""
    with open(path, "rb") as f:
        content = f.read()
    buf = BytesIO(content)
    buf.seek(0)
    return buf
