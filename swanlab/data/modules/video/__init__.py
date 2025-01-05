from swankit.core import MediaBuffer, DataSuite as D, MediaType
from typing import Union, Any, TYPE_CHECKING, Dict, Optional
from io import BytesIO
import os
import tempfile
import secrets
import string


if TYPE_CHECKING:
    # noinspection PyPackageRequirements
    import numpy as _numpy  # type: ignore
    from typing import TextIO

    VideoDataOrPathType = Union["_numpy.ndarray", str, "TextIO", "BytesIO"]
    

def generate_video_filename(format: str) -> str:
    # 生成临时文件路径
    MEDIA_TMP = tempfile.TemporaryDirectory("swanlab-media")
    random_id = generate_id()
    filename = os.path.join(
        MEDIA_TMP.name, random_id + "." + format
    )
    return filename

def generate_id(length: int = 8) -> str:
    """生成一个随机base-36字符串"""
    # There are ~2.8T base-36 8-digit strings. If we generate 210k ids,
    # we'll have a ~1% chance of collision.
    alphabet = string.ascii_lowercase + string.digits
    return "".join(secrets.choice(alphabet) for _ in range(length))


def is_numpy_array(obj: Any) -> bool:
    """检查传入对象是否为numpy数组"""
    return _numpy and isinstance(obj, _numpy.ndarray)


def write_gif_with_image_io(clip: Any, filename: str, fps: int = None) -> None:
    """
    使用imageio库将视频编码为gif格式
    """
    try:
        # noinspection PyPackageRequirements
        import imageio
    except ImportError:
        raise ImportError(
            "swanlab.Video requires imageio when passing raw data. Install with `pip install imageio`"
        )

    writer = imageio.save(filename, fps=clip.fps, quantizer=0, palettesize=256, loop=0)

    for frame in clip.iter_frames(fps=fps, dtype="uint8"):
        writer.append_data(frame)

    writer.close()


class Video(MediaType):
    ACCEPT_FORMAT = ["gif", "mp4", "webm", "ogg"]
    """Video class constructon.
    
    Parameters
    ----------
    data_or_path: (str or BytesIO or numpy.ndarray or torch.Tensor)
        Video can be initialized with a path to a file or an io object.
        The format must be "gif", "mp4", "webm" or "ogg".
        The format must be specified with the format argument.
        Video can be initialized with a numpy tensor.
        The numpy tensor must be either 4 dimensional or 5 dimensional.
        Channels should be (time, channel, height, width) or
        (batch, time, channel, height width)
    caption: (str)
        Caption for the video. Used for display in the SwanLab Dashboard.
    fps: (int)
        The frame rate to use when encoding raw video frames. Default value is 4.
        This parameter has no effect when data_or_path is a string, or bytes.
    format: (str)
        format of video, necessary if initializing with path or io object.
    """
    
    def __init__(
        self,
        data_or_path: "VideoDataOrPathType",
        caption: str = None,
        fps: int = None,
        format: str = None,
    ):
        super().__init__()

        # 检查格式是否正确, 默认格式为gif
        self._format = self.__convert_file_type(format)
        self._width = None
        self._height = None
        self._channels = None
        
        self.buffer = MediaBuffer()

        # 检查fps是否正确, 当提供文件路径或原始字节时，fps参数不影响视频的帧率
        if isinstance(data_or_path, (BytesIO, str)) and fps:
            msg = (
                "`fps` argument does not affect the frame rate of the video " 
                "when providing a file path or raw bytes."
            )
            print(msg)  # 使用Warning也许更好一些

        # 如果data_or_path是str（文件路径）
        if isinstance(data_or_path, str):
            _, ext = os.path.splitext(data_or_path)
            ext = ext[1:].lower()
            if ext not in Video.ACCEPT_FORMAT:
                raise ValueError("swanlab.Video accepts {} formats".format(", ".join(Video.ACCEPT_FORMAT)))
            # 读取这个文件，写入到BytesIO中
            with open(data_or_path, "rb") as f:
                self.buffer.write(f.read())
        # 如果data_or_path是BytesIO
        elif isinstance(data_or_path, BytesIO):
            filename = generate_video_filename(self._format)
            with open(filename, "wb") as f:
                f.write(data_or_path.read())
        else:
            # 如果data_or_path是numpy或torch tensor
            if hasattr(data_or_path, "numpy"):  # TF data eager tensors
                self.video_data = data_or_path.numpy()
            elif is_numpy_array(data_or_path):
                self.video_data = data_or_path
            # 如果都不是，则报错
            else:
                raise ValueError(
                    "swanlab.Video accepts a file path or numpy like data as input"
                )
            
            # 设置fps
            fps = fps or 4
            self.encode(fps=fps)

        self.caption = D.check_caption(caption)

    def __convert_file_type(self, file_type: str = None):
        """转换file_type，并检测file_type是否正确"""
        file_type = file_type or "gif"

        if file_type not in self.ACCEPT_FORMAT:
            raise ValueError(f"swanlab.Video accepts {', '.join(self.ACCEPT_FORMAT)} formats")

        return file_type
    

    def encode(self, fps: int = 4) -> None:
        """
        将numpy数组编码为视频文件

        Args:
            fps (int): 视频的帧率，默认为4帧/秒

        处理流程:
        1. 导入moviepy模块用于视频处理
        2. 将numpy数组转换为合适的视频帧格式
        3. 获取视频的高度、宽度和通道数
        4. 使用moviepy创建视频片段，设置fps
        5. 生成临时文件路径保存视频
        6. 根据moviepy版本尝试不同参数写入视频:
           - 先尝试使用logger=None
           - 如果不支持，尝试verbose=False和progress_bar=False
           - 如果仍不支持，仅使用verbose=False
        7. gif格式使用特殊的write_gif_with_image_io()函数
        8. 最后保存生成的视频文件路径
        """
        try:
            # noinspection PyPackageRequirements
            import moviepy.editor as mpy
        except ImportError:
            raise ImportError(
                "swanlab.Video requires moviepy when passing raw data.  Install with `pip install moviepy`"
            )
            
        # 准备视频数据
        tensor = self._prepare_video(self.video_data)
        _, self._height, self._width, self._channels = tensor.shape  # type: ignore

        # 将图像序列编码为视频
        clip = mpy.ImageSequenceClip(list(tensor), fps=fps)

        # 生成临时文件路径
        filename = generate_video_filename(self._format)
        
        if TYPE_CHECKING:
            kwargs: Dict[str, Optional[bool]] = {}
            
        # 尝试不同的moviepy参数组合写入视频
        try:  # 新版本moviepy支持logger参数
            kwargs = {"logger": None}
            if self._format == "gif":
                write_gif_with_image_io(clip, filename)
            else:
                clip.write_videofile(filename, **kwargs)
        except TypeError:
            try:  # 较老版本moviepy支持progress_bar参数
                kwargs = {"verbose": False, "progress_bar": False}
                if self._format == "gif":
                    clip.write_gif(filename, **kwargs)
                else:
                    clip.write_videofile(filename, **kwargs)
            except TypeError:  # 最老版本moviepy只支持verbose参数
                kwargs = {
                    "verbose": False,
                }
                if self._format == "gif":
                    clip.write_gif(filename, **kwargs)
                else:
                    clip.write_videofile(filename, **kwargs)
                    
        # 将生成的视频文件写入buffer
        with open(filename, "rb") as f:
            self.buffer.write(f.read())

    def _prepare_video(self, video: "_numpy.ndarray") -> "_numpy.ndarray":
        """This logic was mostly taken from tensorboardX."""
        """
        将numpy数组转换为合适的视频帧格式
        
        处理流程:
        1. 加载numpy模块
        2. 检查视频的维度是否至少为4
        3. 如果视频维度小于4，则抛出错误
        4. 如果视频维度为4，则将视频的维度调整为1, *video.shape
        5. 获取视频的batch大小、时间步数、通道数、高度和宽度
        6. 如果视频的数据类型不是uint8，则将视频数据类型转换为uint8
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "swanlab.Video requires numpy when passing raw data. To get it, run `pip install numpy`"
            )
        
        if video.ndim < 4:
            raise ValueError(
                "Video must be atleast 4 dimensions: time, channels, height, width"
            )
        if video.ndim == 4:
            video = video.reshape(1, *video.shape)
        b, t, c, h, w = video.shape

        if video.dtype != np.uint8:
            # logging.warning("Converting video data to uint8")
            video = video.astype(np.uint8)

        def is_power2(num: int) -> bool:
            return num != 0 and ((num & (num - 1)) == 0)

        # 一次性填充到最近的2的幂次方
        if not is_power2(video.shape[0]):
            len_addition = int(2 ** video.shape[0].bit_length() - video.shape[0])
            video = np.concatenate(
                (video, np.zeros(shape=(len_addition, t, c, h, w))), axis=0
            )

        n_rows = 2 ** ((b.bit_length() - 1) // 2)
        n_cols = video.shape[0] // n_rows

        video = video.reshape(n_rows, n_cols, t, c, h, w)
        video = np.transpose(video, axes=(2, 0, 4, 1, 5, 3))
        video = video.reshape(t, n_rows * h, n_cols * w, c)
        return video

    # ---------------------------------- 覆写方法 ----------------------------------

    def parse(self):
        # 文件名称
        hash_name = D.get_hash_by_bytes(self.buffer.getvalue())
        save_name = f"video-step{self.step}-{hash_name}.{self.format}"
        return save_name, self.buffer

    def get_more(self):
        """返回more数据"""
        return {"caption": self.caption} if self.caption is not None else None

    def get_section(self) -> str:
        """设定分组名"""
        return "Video"

    def get_chart(self):
        return self.Chart.VIDEO
