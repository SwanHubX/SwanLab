from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Type, Union

from swanlab.toolkit import MediaType
from .model3d import Model3D
from .molecule import Molecule
from .point_cloud import Box, PointCloud

try:
    # noinspection PyPackageRequirements
    import numpy as np

    ndarray = np.ndarray
except ImportError:
    np = None

    ndarray = None

try:
    from rdkit.Chem import Mol
except ImportError:
    Mol = None


class Object3D:
    """A dispatcher class that converts different types of 3D data to MediaType objects.

    :doc: `Object3D API Reference <https://docs.swanlab.cn/api/py-object3d.html>`

    This class provides a unified interface for handling various 3D data formats including:
    - Point clouds from numpy arrays
    - Point clouds with boxes from {"points": ndarray, "boxes": List[Box]}
    - 3D models from files (.glb)
    - Point cloud files (.swanlab.pts.json)

    Usage Examples:
    --------------
    1. Creating from numpy array (point cloud):
        >>> import numpy as np
        >>> # XYZ coordinates (N, 3)
        >>> points_xyz = np.random.rand(100, 3)
        >>> obj1 = Object3D(points_xyz)
        >>>
        >>> # XYZC coordinates with category (N, 4)
        >>> points_xyzc = np.random.rand(100, 4)
        >>> obj2 = Object3D(points_xyzc)
        >>>
        >>> # XYZRGB coordinates with colors (N, 6)
        >>> points_xyzrgb = np.random.rand(100, 6)
        >>> obj3 = Object3D(points_xyzrgb, caption="Colored Points")

    2. Creating from files:
        >>> # From GLB file
        >>> obj4 = Object3D("model.glb")
        >>>
        >>> # From SwanLab point cloud file
        >>> obj5 = Object3D("points.swanlab.pts.json")

    3. With optional parameters:
        >>> obj6 = Object3D(
        ...     points_xyz,
        ...     caption="My 3D"   # Caption text
        ... )

    4. Creating PointCloud With boxes:
        >>> obj7 = Object3D(
        ...     {"points": points_xyz, "boxes": list(Box)}
        ... )

    5. Creating from Molecule:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> molecule = Molecule.from_mol(mol, caption="Ethanol")
        >>> obj8 = Object3D(molecule)

    Args:
        data: Input data, can be:
            - numpy.ndarray: Point cloud data with shape (N, C) where C is 3,4 or 6
            - str/Path: Path to a 3D file (.glb or .swanlab.pts.json)
            - Molecule: A Molecule object
        caption: Optional description text
        **kwargs: Additional keyword arguments passed to specific handlers

    Returns:
        MediaType: A MediaType object (PointCloud or Model3D)

    Raises:
        ImportError: If numpy is not installed
        ValueError: If input format is not supported
        FileNotFoundError: If input file does not exist
        TypeError: If input type is not supported
    """

    def __new__(
        cls,
        data: Union[ndarray, str, Path, Dict, Mol],
        *,
        caption: Optional[str] = None,
        **kwargs,
    ) -> MediaType:
        cls._check_numpy()

        kwargs['caption'] = caption

        if isinstance(data, ndarray):
            return cls._handle_ndarray(data, **kwargs)

        if isinstance(data, (str, Path)):
            return cls._handle_file(Path(data), **kwargs)

        if isinstance(data, Mol):
            return Molecule.from_mol(data, **kwargs)

        return cls._handle_data(data, **kwargs)

    @staticmethod
    def _check_numpy() -> None:
        if np is None:
            raise ImportError("Numpy is required for Object3D class. Please install it with: pip install numpy.")

    @classmethod
    def _handle_ndarray(cls, data: ndarray, **kwargs) -> MediaType:
        point_cloud_channel_handlers = {
            3: [PointCloud.from_xyz],  # xyz
            4: [PointCloud.from_xyzc],  # xyzc
            6: [PointCloud.from_xyzrgb],  # xyzrgb
        }

        if data.ndim == 2 and data.shape[1] in point_cloud_channel_handlers:
            return cls._try_all(point_cloud_channel_handlers[data.shape[1]], data, **kwargs)

        raise ValueError(
            f"Unsupported array format: shape={data.shape}. "
            f"Expected 2D array with {list(point_cloud_channel_handlers.keys())} channels"
        )

    _FILE_HANDLERS: Dict[str, List[Callable]] = {
        '.swanlab.pts.json': [PointCloud.from_swanlab_pts_json_file],
        '.glb': [Model3D.from_glb_file],
        '.sd': [Molecule.from_sdf_file],
        '.sdf': [Molecule.from_sdf_file],
        '.mol': [Molecule.from_mol_file],
        '.pdb': [Molecule.from_pdb_file],
    }

    @classmethod
    def _handle_file(cls, path: Path, **kwargs) -> MediaType:
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
                    return cls._try_all(handler, path, **kwargs)
                except Exception as e:
                    raise ValueError(
                        f"Error processing file {path} with handler for {suffix}:\n{str(e)}",
                    ) from e

        raise ValueError(
            f"Unsupported file type: {path.name}\n"
            f"File extension: {''.join(suffixes)}\n"
            f"Tried extensions: {', '.join(all_tried_suffixes)}\n"
            f"Supported extensions: {', '.join(sorted(cls._FILE_HANDLERS.keys()))}"
        )

    _TYPE_HANDLERS: Dict[Type, List[Callable]] = {
        dict: [PointCloud.from_swanlab_pts],
    }

    @classmethod
    def _handle_data(cls, data: Any, **kwargs) -> MediaType:
        handler = cls._TYPE_HANDLERS.get(type(data))

        if handler is not None:
            return cls._try_all(handler, data, **kwargs)

        raise TypeError(f"Unsupported input type: {type(data)}")

    @classmethod
    def _try_all(cls, handlers: List[Callable], data: Any, **kwargs) -> MediaType:
        errors = []

        for handler in handlers:
            try:
                return handler(data, **kwargs)
            except Exception as e:
                errors.append(f"Error processing data with handler {handler.__name__}: {str(e)}")

        error_message = ';\n'.join(errors)
        raise ValueError(f"All handlers failed:\n{error_message}")

    @classmethod
    def from_point_data(
        cls,
        points: ndarray,
        *,
        boxes: Optional[List[Box]] = None,
        caption: Optional[str] = None,
        **kwargs,
    ):
        return PointCloud.from_swanlab_pts(
            {
                'points': points,
                'boxes': boxes,
            },
            caption=caption,
            **kwargs,
        )
