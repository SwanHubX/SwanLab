from pathlib import Path
from typing import Any, Callable, Dict, Type, Union

from swankit.core.data import MediaType

from .point_cloud import PointCloud

try:
    import numpy as np
except ImportError:
    np = None


class Object3D:
    def __new__(cls, data: Union[np.ndarray, str, Path]) -> MediaType:
        cls._check_numpy()

        if isinstance(data, np.ndarray):
            return cls._handle_ndarray(data)

        if isinstance(data, (str, Path)):
            return cls._handle_file(Path(data))

        return cls._handle_data(data)

    @classmethod
    def _check_numpy(cls) -> None:
        if np is None:
            raise ImportError("Numpy is required for Object3D class. " "Please install it with: pip install numpy.")

    @classmethod
    def _handle_ndarray(cls, data: np.ndarray) -> MediaType:
        point_cloud_channel_handlers = {
            3: PointCloud.from_xyz,  # xyz
            4: PointCloud.from_xyzc,  # xyzc
            6: PointCloud.from_xyzrgb,  # xyzrgb
        }

        if data.ndim == 2 and data.shape[1] in point_cloud_channel_handlers:
            return point_cloud_channel_handlers[data.shape[1]](data)

        raise ValueError(
            f"Unsupported array format: shape={data.shape}. "
            f"Expected 2D array with {list(point_cloud_channel_handlers.keys())} channels"
        )

    _FILE_HANDLERS: Dict[str, Callable] = {
        '.txt': lambda p: PointCloud.from_file(p),
        '.xyz': lambda p: PointCloud.from_file(p),
        '.pts': lambda p: PointCloud.from_file(p),
        '.ply': lambda p: PointCloud.from_file(p),
        '.pcd': lambda p: PointCloud.from_file(p),
        '.pts.json': lambda p: PointCloud.from_file(p),
    }

    @classmethod
    def _handle_file(cls, path: Path) -> MediaType:
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")

        if not path.is_file():
            raise ValueError(f"Path is not a file: {path}")

        suffixes = path.suffixes  # 例如: ['.xyz', '.pts', '.json']
        if not suffixes:
            raise ValueError(f"File has no extension name: {path}")

        all_tried_suffixes = [''.join(suffixes[i:]).lower() for i in range(len(suffixes))]

        for suffix in all_tried_suffixes:
            handler = cls._FILE_HANDLERS.get(suffix)
            if handler is not None:
                try:
                    return handler(path)
                except Exception as e:
                    raise ValueError(f"Error processing file {path} with handler for {suffix}: {str(e)}") from e

        raise ValueError(
            f"Unsupported file type: {path.name}\n"
            f"File extension: {''.join(suffixes)}\n"
            f"Tried extensions: {', '.join(all_tried_suffixes)}\n"
            f"Supported extensions: {', '.join(sorted(cls._FILE_HANDLERS.keys()))}"
        )

    # 数据类型到处理函数的映射
    _TYPE_HANDLERS: Dict[Type, Callable] = {
        # 可以继续添加其他数据类型的处理器
    }

    @classmethod
    def _handle_data(cls, data: Any) -> MediaType:
        handler = cls._TYPE_HANDLERS.get(type(data))

        if handler is not None:
            return handler(data)

        raise TypeError(f"Unsupported input type: {type(data)}")
