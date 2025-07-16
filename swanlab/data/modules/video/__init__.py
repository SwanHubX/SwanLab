from swanlab.data.modules.image import Image
class Video(Image):    
    def __init__(
        self,
        data_or_path: str,
        caption: str = None,
    ):
        super().__init__(data_or_path, caption=caption)