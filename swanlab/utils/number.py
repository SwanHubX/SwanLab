#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------------------------------- 科学计数法 ----------------------------------


# 如果需求简单，只需函数即可完成系统功能，可将该类删除
class ScientificNotation:
    # 转化阈值 —— 最小值
    __min_threshold = 1e-21

    # 转化阈值 —— 最大值
    __max_threshold = 1e10

    def __init__(self, min_threshold: float = 1e-21, max_threshold: float = 1e10):
        self.__min_threshold = min_threshold
        self.__max_threshold = max_threshold

    def format(self, value, digits=6) -> float:
        """将超出阈值的数字转成科学计数法
        为了便于存储，在返回时使用字符串格式

        Parameters
        ----------
        value : float
            需要判断和转化的格式
        digits : int
            科学计数法时的保留位数

        Returns
        -------
        float
            转化后的数字
        """
        if self.__min_threshold <= value <= self.__max_threshold:
            return value
        else:
            return format(value, f".{digits}e")


# 最大最小的转化阈值
__min_number = 1e-7
__max_number = 1e13
# 正常情况下，小数部分的精确位数，因为e-7时会转成科学计数法，所以最多七位小数
__digits = ".7f"


def sl_format(number) -> [float, int]:
    """按照一定规则转化数据
    1. 未超出阈值，返回浮点类型的数据
    2. 超出阈值，返回科学计数法格式的数据
    3. 删除无意义的0

    Parameters
    ----------
    number : [float, int]
        需要判断转化的数据

    Returns
    -------
    [float, int]
        转化后的数据
    """
    if not (isinstance(number, float) or isinstance(number, int)):
        return ValueError("Param type error cause only float and int invalid")

    if __min_number <= number < __max_number:
        number = format(number, __digits)
        stripped_number = number.rstrip("0").rstrip(".")
        return stripped_number
    else:
        parts = format(number, ".4e").split("e")
        base = "{:.4f}".format(float(parts[0]))
        stripped_base = base.rstrip("0").rstrip(".")
        result = stripped_base + "e" + parts[1]
        return result


if __name__ == "__main__":
    num1 = "1.23e-10"
    print(sl_format(num1))

    num2 = 0.00003
    print(sl_format(num2))

    num3 = 0.00000000003
    print(sl_format(num3))

    num4 = 0.000000000000032
    print(sl_format(num4))

    num5 = 1231231231233123
    print(sl_format(num5))

    num6 = 31231233123
    print(sl_format(num6))
