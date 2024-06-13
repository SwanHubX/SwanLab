#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/11 22:06
@File: pytest_watch.py
@IDE: pycharm
@Description:
    测试cli的watch命令
"""
import os
from swanlab.cli import cli
from click.testing import CliRunner
import multiprocessing
import requests
import swanlab
from tutils import TEMP_PATH
from swanlab.env import get_swanlog_dir, reset_env


def mock_swanlog(path=None):
    if path is None:
        path = get_swanlog_dir()
    reset_env()
    swanlab.init(logdir=path, mode="local")
    swanlab.log({"test": 1})
    swanlab.finish()


# 运行任务
# noinspection PyTypeChecker
def runner_watch(args=None):
    if args is None:
        args = []
    runner = CliRunner()
    runner.invoke(cli, ["watch", *args])


# 测试能否ping通
def ping(host="127.0.0.1", port=5092):
    url = f"http://{host}:{port}"  # noqa
    response = requests.get(url, timeout=10)
    assert response.status_code == 200


def test_watch_default():
    """
    测试watch命令的默认情况
    由于全局解释器锁的存在，无法在同一个进程中同时运行两个线程，所以这里使用多进程进行测试
    """

    # 运行任务之前需保证已经存在swanlog文件夹
    mock_swanlog()

    # 运行任务
    p1 = multiprocessing.Process(target=runner_watch)
    p2 = multiprocessing.Process(target=ping)
    p1.start()
    p2.start()
    p2.join()
    p1.kill()


def test_watch_logdir():
    """
    测试watch命令的logdir参数
    由于全局解释器锁的存在，无法在同一个进程中同时运行两个线程，所以这里使用多进程进行测试
    """
    # 新建一个logdir
    logdir = os.path.join(TEMP_PATH, "watch", "swanlog")
    os.makedirs(logdir, exist_ok=True)
    # 运行任务之前需保证已经存在swanlog文件夹
    mock_swanlog(logdir)

    # 运行任务
    p1 = multiprocessing.Process(target=runner_watch, args=["--logdir", logdir])
    p2 = multiprocessing.Process(target=ping)
    p1.start()
    p2.start()
    p2.join()
    p1.kill()
