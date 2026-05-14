"""TensorBoard TFEvent parsing utilities for the SwanLab converter."""

import io
import os
from typing import Any, Dict, List, Tuple

from swanlab import vendor


def find_tfevents(logdir: str, depth: int = 3) -> Dict[str, List[str]]:
    """Recursively find all TFEvent files under *logdir* up to *depth* levels deep.

    :param logdir: Root directory to search.
    :param depth: Maximum directory depth relative to *logdir*.
    :returns: A dict mapping relative directory names to lists of absolute file paths.
    """
    tfevents_dict: Dict[str, List[str]] = {}

    def _get_depth(path: str) -> int:
        return path.count(os.sep)

    base_depth = _get_depth(os.path.abspath(logdir))

    for root, _dirs, files in os.walk(logdir):
        current_depth = _get_depth(os.path.abspath(root)) - base_depth
        if current_depth > depth:
            continue

        for file in files:
            if "tfevents" in file:
                rel_dir = os.path.relpath(root, logdir)
                if rel_dir not in tfevents_dict:
                    tfevents_dict[rel_dir] = []
                tfevents_dict[rel_dir].append(os.path.join(root, file))

    return tfevents_dict


def get_tf_events_tags_type(tf_event_path: str) -> Dict[str, str]:
    """Inspect a single TFEvent file and map each tag to its data type.

    Supported types: ``scalar``, ``image``, ``audio``, ``text``.

    :param tf_event_path: Absolute path to a TFEvent file.
    :returns: Mapping ``{tag: type_str}``.
    """
    tb = vendor.tensorboard
    event_file_loader = tb.backend.event_processing.event_file_loader  # type: ignore
    tensor_util = tb.util.tensor_util  # type: ignore

    tags: Dict[str, str] = {}
    loader = event_file_loader.EventFileLoader(tf_event_path)

    for event in loader.Load():
        for value in event.summary.value:
            if value.tag in tags:
                continue

            # Determine the type based on which protobuf field is populated.
            if value.HasField("simple_value"):
                tags[value.tag] = "scalar"
            elif value.HasField("image"):
                tags[value.tag] = "image"
            elif value.HasField("audio"):
                tags[value.tag] = "audio"
            elif value.HasField("tensor"):
                arr = tensor_util.make_ndarray(value.tensor)
                if arr.dtype.kind in ("U", "S", "O"):
                    tags[value.tag] = "text"
                elif arr.dtype.kind in ("f", "i", "u", "b"):
                    tags[value.tag] = "scalar"

    return tags


def get_tf_events_tags_data(tf_event_path: str, tags: Dict[str, str]) -> Dict[str, List[Tuple[int, Any, int]]]:
    """Extract logged data from a TFEvent file for the given *tags*.

    :param tf_event_path: Absolute path to a TFEvent file.
    :param tags: Mapping ``{tag: type_str}`` produced by :func:`get_tf_events_tags_type`.
    :returns: Mapping ``{tag: [(step, value, wall_time), ...]}``.
              *value* is a raw scalar, a ``PIL.Image.Image``, ``[audio_np, sample_rate]``,
              or a ``str`` depending on the tag type.
    """
    tb = vendor.tensorboard
    event_file_loader = tb.backend.event_processing.event_file_loader  # type: ignore
    tensor_util = tb.util.tensor_util  # type: ignore

    np = vendor.np
    PILImage = vendor.PIL.Image

    tag_data: Dict[str, List[Tuple[int, Any, int]]] = {tag: [] for tag in tags}
    loader = event_file_loader.EventFileLoader(tf_event_path)

    for event in loader.Load():
        wall_time = int(event.wall_time)
        for value in event.summary.value:
            if value.tag not in tag_data:
                continue

            tag_type = tags[value.tag]

            if tag_type == "scalar":
                if value.HasField("simple_value"):
                    tag_data[value.tag].append((event.step, value.simple_value, wall_time))
                elif value.HasField("tensor"):
                    arr = tensor_util.make_ndarray(value.tensor)
                    tag_data[value.tag].append((event.step, float(arr.item()), wall_time))

            elif tag_type == "image" and value.HasField("image"):
                img_bytes = value.image.encoded_image_string
                image = PILImage.open(io.BytesIO(img_bytes))
                tag_data[value.tag].append((event.step, image, wall_time))

            elif tag_type == "audio" and value.HasField("audio"):
                audio = value.audio
                audio_np = np.frombuffer(audio.encoded_audio_string, dtype=np.int16)
                sample_rate = audio.sample_rate
                tag_data[value.tag].append((event.step, [audio_np, int(sample_rate)], wall_time))

            elif tag_type == "text" and value.HasField("tensor"):
                arr = tensor_util.make_ndarray(value.tensor)
                item = arr.item()
                if isinstance(item, bytes):
                    text = item.decode("utf-8")
                else:
                    text = str(item)
                tag_data[value.tag].append((event.step, text, wall_time))

    return tag_data
