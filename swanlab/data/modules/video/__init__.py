from io import BytesIO
from typing import Tuple, Union, Optional, Dict

from swanlab.data.modules.image import Image
from swanlab.toolkit import MediaType


class Video(MediaType):
    """
    Video module. Currently, only GIF format file paths are supported.

    Args:
        data_or_path: The path to the video file.
        caption: The caption of the video.
    """

    def __init__(self, data_or_path: str, caption: str = None):
        super().__init__()
        if not data_or_path.endswith(".gif"):
            raise ValueError("swanlab.Video only supports gif format file paths")

        self._image = Image(data_or_path, caption=caption, file_type="gif")

    def parse(self) -> Tuple[Union[str, float], Optional[BytesIO]]:
        return self._image.parse()

    def get_section(self) -> str:
        return self._image.get_section()

    def get_more(self) -> Optional[Dict]:
        return self._image.get_more()

    def get_chart(self):
        return self._image.get_chart()
