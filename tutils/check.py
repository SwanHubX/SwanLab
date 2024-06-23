#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/22 19:18
@File: check.py.py
@IDE: pycharm
@Description:
    开发/测试运行前检查，并加载.env文件设置的环境变量
    此文件应该在顶部导入
"""
import os
import subprocess
import sys
from dotenv import load_dotenv

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
load_dotenv(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".env"))

is_pytest_env = 'PYTEST_VERSION' in os.environ
is_skip_test = os.getenv("TEST_CLOUD_SKIP") is not None
is_cloud_dev_env = os.getenv("SWANLAB_API_HOST") is not None and os.getenv("SWANLAB_WEB_HOST") is not None
if not is_cloud_dev_env:
    # 测试环境
    if is_pytest_env and not is_skip_test:
        print("请设置开发云服务环境变量，或者设置环境变量TEST_CLOUD_SKIP以跳过云测试", file=sys.stderr)
        sys.exit(2)
    # 开发环境
    elif not is_pytest_env:
        print("请设置开发云服务环境变量以运行开发测试脚本", file=sys.stderr)
        sys.exit(2)
