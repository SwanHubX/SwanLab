"""
@author: cunyue
@file: bimap.py
@time: 2026/3/15 01:04
@description: 双向映射适配器
"""


class BiMap:
    """双向映射适配器"""

    def __init__(self, mapping: dict):
        self._forward = mapping
        self._reverse = {v: k for k, v in mapping.items()}

    def __getitem__(self, key):
        if key in self._forward:
            return self._forward[key]
        if key in self._reverse:
            return self._reverse[key]
        raise KeyError(f"No mapping found for: {key} when accessing {self.__class__.__name__}")

    def get(self, key, default=None):
        return self._forward.get(key) or self._reverse.get(key, default)
