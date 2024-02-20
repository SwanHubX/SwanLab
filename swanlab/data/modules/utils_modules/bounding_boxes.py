import numpy as np
import numbers
from typing import Mapping, Any


class BoundingBoxes:
    """
    用于处理2D边界框（bounding boxes）数据的类。
    """

    # box_item的格式形如：{"box_data":..., "class_labels":...}
    # box_data存储了所有box的数据(包括位置、类别id、准确率、是否有caption、...)，class_labels存储了所有的标签合集

    def __init__(self, box_item: dict, key: str) -> None:
        # 验证box_item是否符合要求
        self.validate(box_item)

        self._box_data = box_item["box_data"]
        self._key = key

        # 如果box_item中不存在class_labels参数, 则生成默认格式的类型标签
        if "class_labels" not in box_item:
            classes = np.unique(list(box["class_id"] for box in box_item["box_data"])).astype(np.int32).tolist()
            class_labels = {c: "class_" + str(c) for c in classes}
            self._class_labels = class_labels
        else:
            self._class_labels = box_item["class_labels"]

    def validate(self, box_item) -> bool:
        """
        检查box_data是否符合要求
        """

        # 如果box_item中包含class_labels参数，则检查参数是否符合规范
        if "class_labels" in box_item:
            for k, v in list(box_item["class_labels"].items()):
                if (not isinstance(k, numbers.Number)) or (not isinstance(v, str)):
                    raise TypeError("Class labels must be a dictionary of numbers to string")

        boxes = box_item["box_data"]

        # 如果boxes不是list类型，则报错
        if not isinstance(boxes, list):
            raise TypeError("'box_data' must be a list")

        # 对于boxes中的每一个box
        for box in boxes:
            error_str = "Each box in box_data must contain a position with: middle, width, and height or \
                    \nminX, maxX, minY, maxY."

            # 验证position参数（必须参数）
            # 如果box中不包含position参数，则报错
            if "position" not in box:
                raise TypeError(error_str)
            else:
                valid = False
                if (
                    "middle" in box["position"]
                    and len(box["position"]["middle"]) == 2
                    and validate_dict_key_has_number(box["position"], "width")
                    and validate_dict_key_has_number(box["position"], "height")
                ):
                    valid = True

                elif (
                    validate_dict_key_has_number(box["position"], "minX")
                    and validate_dict_key_has_number(box["position"], "maxX")
                    and validate_dict_key_has_number(box["position"], "minY")
                    and validate_dict_key_has_number(box["position"], "maxY")
                ):
                    valid = True

                if not valid:
                    raise TypeError(error_str)

            # 验证score参数（可选参数）
            if ("scores" in box) and not isinstance(box["scores"], dict):
                raise TypeError("Box scores must be a dictionary")
            elif "scores" in box:
                for k, v in list(box["scores"].items()):
                    if not isinstance(k, str):
                        raise TypeError("A score key must be a string")
                    if not isinstance(v, numbers.Number):
                        raise TypeError("A score value must be a number")

            # 验证class_id参数
            if ("class_id" in box) and not isinstance(box["class_id"], int):
                raise TypeError("A box's class_id must be an integer")

            # 验证box_caption参数
            if ("box_caption" in box) and not isinstance(box["box_caption"], str):
                raise TypeError("A box's caption must be a string")

        return True


def validate_dict_key_has_number(dictionary: Mapping, key: Any) -> bool:
    return key in dictionary and isinstance(dictionary[key], numbers.Number)
