#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/1/27 10:00
@File: offline.py
@IDE: pycharm
@Description:
    offline命令，设置SWANLAB_MODE环境变量为offline
"""

import os
import click

from swanlab.env import SwanLabEnv


@click.command()
def offline():
    """Set SwanLab mode to offline.
    
    This command sets the SWANLAB_MODE environment variable to "offline",
    which means SwanLab will save data locally without uploading to the cloud.
    
    Example:
        swanlab offline
        python your_script.py  # Now runs in offline mode
    """
    mode_key = SwanLabEnv.MODE.value
    os.environ[mode_key] = "offline"
    click.echo(f"✅ SwanLab mode set to offline (SWANLAB_MODE={os.environ[mode_key]})")
    click.echo("💡 Your next SwanLab experiment will run in offline mode.")
    click.echo("   Data will be saved locally without uploading to the cloud.") 