import numpy as np
import numbers


class ImageMask:
    def __init__(self, mask_item: dict, key: str) -> None:
        # mask_item含有两个参数: mask_data和class_labels
        # mask_data参数的类型为numpy.ndarray

        # 如果mask_item中不存在class_labels参数, 则生成默认格式的类型标签
        if "class_labels" not in mask_item:
            classes = np.unique(list(box["class_id"] for box in mask_item["box_data"])).astype(np.int32).tolist()
            class_labels = {c: "class_" + str(c) for c in classes}
            self._class_labels = class_labels
        else:
            self._class_labels = mask_item["class_labels"]

        # 验证mask_item是否符合要求
        self.validate(mask_item)
        self._mask_data = mask_item["mask_data"]
        self._key = key

    def validate(self, mask_item) -> bool:
        """
        检查mask_data是否符合要求
        """

        # 如果mask_item中包含class_labels参数，则检查参数是否符合规范
        if "class_labels" in mask_item:
            for k, v in list(mask_item["class_labels"].items()):
                if (not isinstance(k, numbers.Number)) or (not isinstance(v, str)):
                    raise TypeError("Class labels must be a dictionary of numbers to string")

        if "mask_data" not in mask_item:
            raise TypeError(
                'swanlab.Image -> masks Missing key "mask_data": An image mask requires mask data: a 2D array representing the predictions'
            )
        else:
            error_str = "mask_data must be a 2D array"
            shape = mask_item["mask_data"].shape
            if len(shape) != 2:
                raise TypeError(error_str)
            if not ((mask_item["mask_data"] >= 0).all() and (mask_item["mask_data"] <= 255).all()) and issubclass(
                mask_item["mask_data"].dtype.type, np.integer
            ):
                raise TypeError("Mask data must be integers between 0 and 255")

        return True
