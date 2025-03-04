from dataclasses import dataclass
from functools import once
from pathlib import Path
from typing import Dict, Optional, Tuple

from swankit.core.data import DataSuite as D
from swankit.core.data import MediaBuffer, MediaType


@dataclass(frozen=True)
class Model3D(MediaType):
    """3D model data representation

    Attributes:
        glb_path: Path to the GLB file
        step: Optional step number for visualization
        caption: Optional description text
    """

    glb_path: Path
    step: Optional[int] = None
    caption: Optional[str] = None

    def __post_init__(self):
        """Validate input data after initialization"""
        if not self.glb_path.exists():
            raise FileNotFoundError(f"GLB file not found: {self.glb_path}")

        if not self.glb_path.is_file():
            raise ValueError(f"Path is not a file: {self.glb_path}")

        if self.glb_path.suffix.lower() != '.glb':
            raise ValueError(f"File must be a GLB file: {self.glb_path}")

    @classmethod
    def from_glb_file(
        cls, path: Path, *, step: Optional[int] = None, caption: Optional[str] = None, **kwargs
    ) -> "Model3D":
        """Create Model3D from GLB file

        Args:
            path: Path to the .glb file
            step: Optional step number
            caption: Optional description text

        Returns:
            Model3D object

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file is not a GLB file
        """
        return cls(path, step=step, caption=caption, **kwargs)

    # ---------------------------------- override ----------------------------------

    @once
    def parse(self) -> Tuple[str, MediaBuffer]:
        """Convert model to buffer for transmission"""
        # 直接复制GLB文件内容到buffer
        with open(self.glb_path, 'rb') as f:
            buffer = MediaBuffer()
            buffer.write(f.read())

        buffer.seek(0)
        hash_name = D.get_hash_by_bytes(buffer.getvalue())[:16]
        save_name = f"{self.glb_path.stem}-step{self.step}-{hash_name}.glb"
        buffer.seek(0)

        return save_name, buffer

    def get_chart(self) -> MediaType.Chart:
        return MediaType.Chart.OBJECT3D

    def get_section(self) -> str:
        return "Model3d"

    def get_more(self) -> Optional[Dict[str, str]]:
        return {"caption": self.caption} if self.caption else None
