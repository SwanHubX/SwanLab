import os
from PIL import Image
import io
import numpy as np

try:
    import tensorflow as tf
except ImportError as e:
    raise TypeError(
        "Tensorboard Converter requires tensorflow when process tfevents file. Install with 'pip install tensorflow'."
    )


def get_tf_events_tags_type(tf_event_path: str):
    """获取TFEvent文件中所有tag的类型，并返回一个字典
    比如{"tag1": "scalar", "tag2": "image", "tag3": "audio", "tag4": "text"}

    Args:
        tf_event_path (_type_): 单个tf_event文件的路径

    Returns:
        _type_: 返回一个字典，键是tag，值是tag的类型
    """
    # 确保路径存在
    assert os.path.exists(tf_event_path), "TFEvent file does not exist"

    # 用来存储所有tag的集合和类型
    tags = {}

    # 遍历所有事件
    for event in tf.compat.v1.train.summary_iterator(tf_event_path):
        for value in event.summary.value:
            if value.tag not in tags:
                if value.HasField("simple_value"):
                    tags[value.tag] = "scalar"
                elif value.HasField("image"):
                    tags[value.tag] = "image"
                elif value.HasField("audio"):
                    tags[value.tag] = "audio"
                elif value.HasField("tensor") and value.tensor.dtype == tf.string:
                    tags[value.tag] = "text"

    return tags


def get_tf_events_tags_data(tf_event_path: str, tags: dict):
    """获取TFEvent文件中所有tag的数据，并返回一个字典
    比如{"tag1": [(step1, value1), (step2, value2)], "tag2": [(step1, value1), (step2, value2)]}

    Args:
        tf_event_path (str): tf_event文件的路径
        tags (dict): tag的类型字典, 形如{"loss": "scalar", "np_image": "image"}

    Returns:
        _type_: 返回一个字典，键是tag，值是完整的数据, 形如{"loss": [(step1, value1), (step2, value2)], "np_image": [(step1, value1), (step2, value2)]}
    """

    # 用来存储每个tag的数据
    tag_data = {tag: [] for tag in tags}

    # 再次遍历文件，这次是为了提取数据
    for event in tf.compat.v1.train.summary_iterator(tf_event_path):
        # wall_time = datetime.datetime.fromtimestamp(event.wall_time).strftime("%Y-%m-%d %H:%M:%S")
        wall_time = int(event.wall_time)
        for value in event.summary.value:
            if value.tag in tag_data:
                if tags[value.tag] == "scalar":
                    tag_data[value.tag].append((event.step, value.simple_value, wall_time))
                elif tags[value.tag] == "image" and value.HasField("image"):
                    img_str = value.image.encoded_image_string
                    image = Image.open(io.BytesIO(img_str))
                    tag_data[value.tag].append((event.step, image, wall_time))
                elif tags[value.tag] == "audio" and value.HasField("audio"):
                    audio = value.audio
                    audio_np = np.frombuffer(value.audio.encoded_audio_string, dtype=np.int16)
                    sample_rate = audio.sample_rate

                    tag_data[value.tag].append((event.step, [audio_np, int(sample_rate)], wall_time))
                elif tags[value.tag] == "text" and value.HasField("tensor"):
                    text = tf.make_ndarray(value.tensor).item().decode("utf-8")
                    tag_data[value.tag].append((event.step, text, wall_time))

    return tag_data


def find_tfevents(logdir: str, depth: int = 3):
    """查找指定目录下的所有tfevents文件

    Args:
        logdir (str): 日志文件夹路径
        depth (int, optional): 目录深度，默认为3

    Returns:
        Dict: 返回一个字典，键是子目录名，值是该目录下的所有tfevents文件路径列表
    """
    tfevents_dict = {}

    def get_depth(path):
        return path.count(os.sep)

    base_depth = get_depth(logdir)

    for root, dirs, files in os.walk(logdir):
        current_depth = get_depth(root) - base_depth
        if current_depth > depth:
            continue

        for file in files:
            if "tfevents" in file:
                # 使用目录名作为字典的键
                directory_key = os.path.relpath(root, logdir)  # 获取目录的最后一个部分作为键
                if directory_key not in tfevents_dict:
                    tfevents_dict[directory_key] = []
                tfevents_dict[directory_key].append(os.path.join(root, file))

    return tfevents_dict
