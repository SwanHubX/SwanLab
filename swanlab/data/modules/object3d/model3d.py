from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

from swanlab.toolkit import DataSuite as D, MediaBuffer, MediaType


@dataclass()
class Model3D(MediaType):
    """3D model data representation class that handles GLB format files.

    This class provides functionality to handle 3D model data in GLB format for visualization.
    It supports loading GLB files and manages metadata like caption.

    Attributes:
        glb_path: Path to the GLB file
        caption: Optional description text

    Examples:
        >>> # Create from GLB file
        >>> model = Model3D.from_glb_file("model.glb")
        >>>
        >>> # Create with metadata
        >>> model = Model3D.from_glb_file(
        ...     "model.glb",
        ...     caption="My 3D Model"
        ... )
        >>>
        >>> # Use with SwanLab
        >>> import swanlab
        >>> swanlab.log({"model": model})
    """

    glb_path: Path
    step: Optional[int] = None
    caption: Optional[str] = None

    def __post_init__(self):
        """Validate input data after initialization.

        Checks:
            1. File exists
            2. Path points to a file
            3. File has .glb extension

        Raises:
            FileNotFoundError: If GLB file does not exist
            ValueError: If path is not a file or not a GLB file
        """
        if not self.glb_path.exists():
            raise FileNotFoundError(f"GLB file not found: {self.glb_path}")

        if not self.glb_path.is_file():
            raise ValueError(f"Path is not a file: {self.glb_path}")

        if self.glb_path.suffix.lower() != '.glb':
            raise ValueError(f"File must be a GLB file: {self.glb_path}")

    @classmethod
    def from_glb_file(cls, path: Path, *, caption: Optional[str] = None, **kwargs) -> "Model3D":
        """Create Model3D instance from a GLB file.

        This is the main factory method to create Model3D objects. It handles file validation
        and metadata management.

        Args:
            path: Path to the .glb file
            caption: Optional description text for the model
            **kwargs: Additional keyword arguments passed to constructor

        Returns:
            Model3D: A new Model3D instance

        Raises:
            FileNotFoundError: If GLB file does not exist
            ValueError: If file is not a GLB file

        Examples:
            >>> # Basic usage
            >>> model = Model3D.from_glb_file("model.glb")
            >>>
            >>> # With metadata
            >>> model = Model3D.from_glb_file(
            ...     "model.glb",
            ...     caption="My 3D Model"
            ... )
        """
        return cls(path, caption=caption, **kwargs)

    # ---------------------------------- override ----------------------------------

    def parse(self) -> Tuple[str, MediaBuffer]:
        """Convert model to buffer for transmission.

        This method reads the GLB file and prepares it for transmission. It:
        1. Reads the file content into a buffer
        2. Generates a unique hash from the content
        3. Creates a unique filename using the original name, step, and hash

        Returns:
            Tuple[str, MediaBuffer]: (filename, file_content)
            - filename: Generated unique filename
            - file_content: File content in MediaBuffer

        Example:
            >>> model = Model3D.from_glb_file("chair.glb")
            >>> filename, content = model.parse()
            >>> print(filename)
            'chair-step1-a1b2c3d4e5f6g7h8.glb'
        """
        with open(self.glb_path, 'rb') as f:
            buffer = MediaBuffer()
            buffer.write(f.read())

        buffer.seek(0)
        hash_name = D.get_hash_by_bytes(buffer.getvalue())[:16]
        save_name = f"{self.glb_path.stem}-step{self.step}-{hash_name}.glb"
        buffer.seek(0)

        return save_name, buffer

    def get_chart(self) -> MediaType.Chart:
        """Get chart type for visualization.

        Returns:
            MediaType.Chart.OBJECT3D: Indicates this is a 3D object
        """
        return MediaType.Chart.OBJECT3D

    def get_section(self) -> str:
        """Get section name for organization.

        Returns:
            str: Section name "Model3d"
        """
        return "Model3d"

    def get_more(self) -> Optional[Dict[str, str]]:
        """Get additional metadata.

        Returns:
            Optional[Dict[str, str]]: Dictionary containing caption if set
        """
        return {"caption": self.caption} if self.caption else None
