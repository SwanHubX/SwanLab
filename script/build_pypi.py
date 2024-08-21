#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-11 13:09:58
@File: build_pypi.py
@IDE: vscode
@Description:
    构建项目的脚本，用于构建pypi包
"""
import subprocess
import shutil
import json
import os

# 如果dist文件夹存在则删除
if os.path.exists("../dist"):
    shutil.rmtree("../dist")

# 设置版本号
version = os.getenv("VERSION")
if not version:
    raise ValueError("尚未指定构建版本号，请设置VERSION环境变量，如`export VERSION=0.6.0`")
with open("swanlab/package.json", 'r+') as f:
    p = json.load(f)
    p["version"] = version
    f.seek(0)
    json.dump(p, f, indent=4)
    f.truncate()

# 构建python项目
subprocess.run("python -m build", shell=True)
