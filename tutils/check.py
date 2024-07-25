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
swanboard_version = [i.split(": ")[1] for i in swanboard.split("\n") if i.startswith("Version")][0].split("\r")[0]
swankit_version = [i.split(": ")[1] for i in swankit.split("\n") if i.startswith("Version")][0].split("\r")[0]
with open(os.path.join(swanlab_dir, "requirements.txt"), "r") as f:
    packages = f.read().split("\n")
packages = [i for i in packages if "swanboard" in i or "swankit" in i]
for i in packages:
    if "swanboard" in i and swanboard_version not in i:
        raise Exception(f"swanboard过时，运行 pip install -r requirements.txt 进行更新.")
    if "swankit" in i and swankit_version not in i:
        raise Exception(f"swankit过时，运行 pip install -r requirements.txt 进行更新.")

# ---------------------------------- 检查是否跳过云测试 ----------------------------------
load_dotenv(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".env"))
runtime = os.getenv("SWANLAB_RUNTIME")
# pytest测试环境
is_pytest_env = "PYTEST_VERSION" in os.environ
# 是否跳过部分云端测试
is_skip_cloud_test = runtime == 'test-no-cloud'
# 是否为测试环境
is_test_runtime = os.getenv("SWANLAB_RUNTIME") in ['test', 'test-no-cloud']
# 如果为pytest测试环境，环境变量SWANLAB_RUNTIME必须为['test', 'test-no-cloud']之一
# 如果没有跳过部分云端测试，必须设置SWANLAB_WEB_HOST、SWANLAB_API_HOST、SWANLAB_API_KEY
"""
* 推荐在项目根目录下设置.env文件完成环境变量的设置，具体代码为：
    
    SWANLAB_RUNTIME=test-no-cloud

* 如果使用终端，也可以根据不同操作系统版本选择需要的命令

    WINDOWS CMD COMMAND: set SWANLAB_RUNTIME="test-no-cloud"
    WINDOWS POWERSHELL COMMAND: $env:SWANLAB_RUNTIME="test-no-cloud"
    MAC & LINUX COMMAND: export SWANLAB_RUNTIME="test-no-cloud"
"""
if is_pytest_env:
    if not is_test_runtime:
        print("请设置SWANLAB_RUNTIME环境变量为 test 或 test-no-cloud 以运行云测试", file=sys.stderr)
        sys.exit(2)
    if not is_skip_cloud_test:
        envs = ["SWANLAB_WEB_HOST", "SWANLAB_API_HOST", "SWANLAB_API_KEY"]
        if not all([os.getenv(i) for i in envs]):
            print("请设置云测试相关环境变量以运行云测试", file=sys.stderr)
            sys.exit(2)
