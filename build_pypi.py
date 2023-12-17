#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-12-11 13:09:58
@File: build_pypi.py
@IDE: vscode
@Description:
    构建项目
"""
import subprocess
import shutil
import os

# 构建node项目
subprocess.run("npm run build.release", shell=True)
# 如果dist文件夹存在则删除
if os.path.exists("dist"):
    shutil.rmtree("dist")
# 构建python项目
subprocess.run("python -m build", shell=True)
