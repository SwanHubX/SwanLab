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

    def __init__(self, data_or_path: Union[str, "Video"], caption: str = None):
        super().__init__()
        # Support swanlab.Video as input (e.g., swanlab.Video(swanlab.Video(path)))
        if isinstance(data_or_path, Video):
            # Use Image's nested support to create a new Image from the existing one
            self._image = Image(
                data_or_path._image,
                caption=caption if caption is not None else data_or_path._image.caption,
                file_type="gif"
            )
            return

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
