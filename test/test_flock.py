#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-01 22:58:44
@File: test\test_flock.py
@IDE: vscode
@Description:
    测试文件锁，用于多进程写入同一个文件
"""
import os
import portalocker
import time
import random

# 生成随机整数
id = random.randint(0, 1000)


def flock_test():
    print("open file")
    # 具有读写权限
    with open("test.txt", "r+") as f:
        print("try to write")
        portalocker.lock(f, portalocker.LOCK_EX)  # 加锁
        content = f.read()
        print("now content is: " + content)
        print("writing")
        f.write("hello world from process" + str(id) + "\n")
        time.sleep(10)
        f.close()  # 关闭自动解锁
    print("done")
    time.sleep(1)


if __name__ == "__main__":
    print(f"process {id} start")

    while True:
        flock_test()
