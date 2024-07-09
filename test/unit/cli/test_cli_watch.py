#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/11 22:06
@File: pytest_watch.py
@IDE: pycharm
@Description:
    测试cli的watch命令
"""
# import os
# import time
# import pytest
# from swanlab.cli import cli
# from click.testing import CliRunner
# import multiprocessing
# import requests
# import swanlab
# from tutils import TEMP_PATH
# from swanlab.env import get_swanlog_dir, SwanLabEnv
#
#
# def mock_swanlog(path=None):
#     if path is None:
#         path = get_swanlog_dir()
#     os.makedirs(path, exist_ok=True)
#     del os.environ[SwanLabEnv.SWANLOG_FOLDER.value]
#     swanlab.init(logdir=path, mode="local")
#     swanlab.log({"test": 1})
#     swanlab.finish()
#
#
# # 运行任务
# # noinspection PyTypeChecker
# def runner_watch(*args):
#     runner = CliRunner()
#     return runner.invoke(cli, ["watch", *args])
#
#
# def runner_watch_wrapper(*args):
#     from tutils import reset_some_env
#     reset_some_env()
#     r = runner_watch(*args)
#     if r.exit_code != 0:
#         raise Exception(r.output)
#
#
# # 测试能否ping通
# def ping(host="127.0.0.1", port=5092):
#     url = f"http://{host}:{port}"  # noqa
#     time.sleep(3)
#     try:
#         response = requests.get(url, timeout=10)
#         assert response.status_code == 200
#     except Exception as e:
#         raise Exception("ping failed", str(e))
#
#
# @pytest.mark.parametrize("logdir, ping_args, args", [
#     # 无参数
#     [None, [], []],
#     # 指定logdir
#     [os.path.join(TEMP_PATH, "watch", "swanlog"), [], ["--logdir", os.path.join(TEMP_PATH, "watch", "swanlog")]],
#     # 指定host和port
#     [None, ["0.0.0.0", "5093"], ["--host", "0.0.0.0", "--port", "5093"]],
#     # 直接watch
#     [os.path.join(TEMP_PATH, "watch", "swanlog"), [], [os.path.join(TEMP_PATH, "watch", "swanlog")]],
# ])
# def test_watch_ok(logdir, ping_args, args):
#     """
#     测试watch命令，正常情况
#     """
#     mock_swanlog(logdir)
#     p1 = multiprocessing.Process(target=runner_watch_wrapper, args=args)
#     p2 = multiprocessing.Process(target=ping, args=ping_args)
#     p1.start()
#     p2.start()
#     p2.join()
#     p1.kill()
#     assert p2.exitcode == 0
#
#
# def test_watch_wrong_logdir():
#     """
#     测试watch命令，logdir不存在
#     """
#     result = runner_watch("wrong_logdir")
#     assert result.exit_code == 2
#     result = runner_watch("--logdir", "wrong_logdir")
#     assert result.exit_code == 2
#     os.environ[SwanLabEnv.SWANLOG_FOLDER.value] = "wrong_logdir"
#     result = runner_watch()
#     assert result.exit_code == 2
#
#
# def test_watch_wrong_host_port():
#     """
#     测试watch命令，host和port错误
#     """
#     mock_swanlog()
#     result = runner_watch("--host", "")
#     assert result.exit_code == 6
#     result = runner_watch("--port", "0")
#     assert result.exit_code == 2
#     result = runner_watch("--port", "65536")
#     assert result.exit_code == 2
#     # 如果ip被占用，会报错
#     p1 = multiprocessing.Process(target=runner_watch_wrapper, args=[get_swanlog_dir(), "--port", "5092"])
#     p1.start()
#     time.sleep(3)
#     result = runner_watch("--port", "5092")
#     assert result.exit_code == 7
#     p1.kill()
