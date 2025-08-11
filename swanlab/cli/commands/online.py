#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/1/27 10:00
@File: online.py
@IDE: pycharm
@Description:
    online命令，在swanlog文件夹中写入设置文件，设置mode为cloud
"""

import os
import json
import click

from swanlab.env import SwanLabEnv


@click.command()
def online():
    """Set SwanLab mode to online (cloud).
    
    This command creates a settings file in the swanlog directory with mode=cloud,
    which means SwanLab will save data locally and upload to the cloud.
    
    Example:
        swanlab online
        python your_script.py  # Now runs in online mode
    """
    # 获取swanlog目录路径
    env_key = SwanLabEnv.SWANLOG_FOLDER.value
    logdir = os.environ.get(env_key) or os.path.join(os.getcwd(), "swanlog")
    logdir = os.path.abspath(logdir)
    
    # 确保swanlog目录存在
    try:
        os.makedirs(logdir, exist_ok=True)
        if not os.access(logdir, os.W_OK):
            raise IOError(f"no write permission for path: {logdir}")
    except Exception as error:
        raise IOError(f"Failed to create or access logdir: {logdir}, error: {error}")
    
    # 创建设置文件路径
    settings_file = os.path.join(logdir, ".swanlab_settings.json")
    
    # 写入设置文件
    settings = {"mode": "cloud"}
    try:
        with open(settings_file, "w", encoding="utf-8") as f:
            json.dump(settings, f, indent=2, ensure_ascii=False)
    except Exception as error:
        raise IOError(f"Failed to write settings file: {settings_file}, error: {error}")
    
    click.echo(f"✅ SwanLab mode set to online (cloud)")
    click.echo(f"📁 Settings file created: {settings_file}")
    click.echo("💡 Your next SwanLab experiment will run in online mode.")
    click.echo("   Data will be saved locally and uploaded to the cloud.") 