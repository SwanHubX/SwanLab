"""
@author: cunyue
@file: main.py
@time: 2026/3/5 13:11
@description: hello world
"""

from typing import TypedDict


def main():
    t: Test = {"a": 2, "b": 2}
    print("hello world", t["a"], t["b"])


class Test(TypedDict):
    a: int
    b: int


if __name__ == "__main__":
    main()
