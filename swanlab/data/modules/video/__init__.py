from swanlab.data.modules.image import Image

class Video(Image):    
    """
    Video module. Currently, only GIF format file paths are supported.

    Args:
        data_or_path: The path to the video file.
        caption: The caption of the video.
    """
    def __init__(
        self,
        data_or_path: str,
        caption: str = None,
    ):
        if not data_or_path.endswith(".gif"):
            raise ValueError("swanlab.Video only supports gif format")
        
        super().__init__(data_or_path, caption=caption, file_type="gif")
    
