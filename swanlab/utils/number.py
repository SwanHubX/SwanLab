#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------------------------------- 科学计数法 ----------------------------------


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


if __name__ == "__main__":
    num1 = "1.23e-10"
    print(num1)
    print(format(float(num1), ".10f"))

    num2 = 123321312312.3123123123123
    print(format(num2, ".12e"))
