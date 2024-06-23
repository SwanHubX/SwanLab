#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/22 19:18
@File: check.py.py
@IDE: pycharm
@Description:
    用于检查requirements.txt中的包是否都已经安装
"""
import os
import subprocess
import sys

swanlab_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------- 检查swanboard、swankit包的版本号与当前系统是否一致 ----------------------------------

swanboard = subprocess.run("pip show swanboard", shell=True, capture_output=True).stdout.decode()
swankit = subprocess.run("pip show swankit", shell=True, capture_output=True).stdout.decode()
swanboard_version = [i.split(": ")[1] for i in swanboard.split("\n") if i.startswith("Version")][0]
swankit_version = [i.split(": ")[1] for i in swankit.split("\n") if i.startswith("Version")][0]
with open(os.path.join(swanlab_dir, "requirements.txt"), "r") as f:
    packages = f.read().split("\n")
packages = [i for i in packages if "swanboard" in i or "swankit" in i]
for i in packages:
    if "swanboard" in i and swanboard_version not in i:
        print(f'swanboard过时，运行 pip install -r requirements.txt 进行更新.', file=sys.stderr)
        sys.exit(2)
    if "swankit" in i and swankit_version not in i:
        print(f'swankit过时，运行 pip install -r requirements.txt 进行更新.', file=sys.stderr)
        sys.exit(2)

# ---------------------------------- 检查是否跳过云测试，如果没跳过，相关环境变量需要指定----------------------------------
is_pytest_env = 'PYTEST_CURRENT_TEST' in os.environ
is_skip_test = os.getenv("TEST_CLOUD_SKIP") is not None
