from swanlab.data.modules.image import Image
from typing import Union

class Video(Image):    
    def __init__(
        self,
        data_or_path: str,
        caption: str = None,
        mode: str = None,
        file_type: str = None,
        size: Union[int, list, tuple] = None,
    ):
        super().__init__(data_or_path, caption=caption, mode=mode, file_type=file_type, size=size)