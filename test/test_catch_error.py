#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-16 14:35:04
@File: test/test_catch_error.py
@IDE: vscode
@Description:
    测试异常抛出
"""
# excepthook.py
import sys, traceback
from datetime import datetime
import atexit

fError = open("except_error.log", "a")


def UserExceptHook(tp, val, tb):
    traceList = traceback.format_tb(tb)
    html = repr(tp) + "\n"
    html += repr(val) + "\n"
    for line in traceList:
        html += line + "\n"
    print(html, file=sys.stderr)
    print(datetime.now(), file=fError)
    print(html, file=fError)
    fError.close()
    import time

    print("执行异常回调函数")
    time.sleep(10)


def close():
    print("执行程序结束的回调.")


def main():
    sFirst = input("First number:")
    sSecond = input("Second number:")
    try:
        fResult = int(sFirst) / int(sSecond)
    except Exception:
        print("发现异常，但我不处理，抛出去.")
        raise
    else:
        print(sFirst, "/", sSecond, "=", fResult)


atexit.register(close)

sys.excepthook = UserExceptHook
main()
fError.close()
