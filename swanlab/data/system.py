#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2023-11-26 16:49:47
@File: swanlab\database\system.py
@IDE: vscode
@Description:
    采集系统数据，包括内存、CPU、GPU、硬盘、网络等
"""
import platform
import socket
import sys
import subprocess
import multiprocessing
import pynvml


def __replace_second_colon(input_string, replacement):
    """替换字符串中第二个‘:’"""
    first_colon_index = input_string.find(":")
    if first_colon_index != -1:
        second_colon_index = input_string.find(":", first_colon_index + 1)
        if second_colon_index != -1:
            return input_string[:second_colon_index] + replacement + input_string[second_colon_index + 1 :]
    return input_string


def __get_remote_url():
    """获取git仓库的远程链接

    Returns
    -------
    str
        git remote url
    """
    try:
        # 运行git命令获取远程仓库URL
        result = subprocess.run(
            ["git", "config", "--get", "remote.origin.url"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        # 检查命令是否成功运行
        if result.returncode == 0:
            url = result.stdout.strip().replace("git@", "https://")
            if url.endswith(".git"):
                url = url[:-4]
            return __replace_second_colon(url, "/")
        else:
            print(f"Error: {result.stderr}")
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def __get_git_branch_and_commit():
    """获取项目git的分支名和该分支下最新提交的hash"""
    try:
        # 获取当前分支名称
        branch_process = subprocess.Popen(
            ["git", "branch", "--show-current"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        branch_name, branch_err = branch_process.communicate()
        if branch_process.returncode != 0:
            print("Error getting branch name:", branch_err)
            return None, None

        branch_name = branch_name.strip().decode("utf-8")

        # 获取当前分支的最新提交hash
        commit_process = subprocess.Popen(
            ["git", "rev-parse", branch_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        commit_hash, commit_err = commit_process.communicate()

        # 如果无法获取最新提交hash, 那么就返回branch_name和None
        if commit_process.returncode != 0:
            print("Error getting commit hash:", commit_err)
            return branch_name, None

        commit_hash = commit_hash.strip().decode("utf-8")
        return branch_name, commit_hash

    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None


def __get_gpu_info():
    """获取 GPU 信息"""
    info = {"cores": None, "type": []}
    try:
        pynvml.nvmlInit()
    except:
        return info
    try:
        # 获取 NVIDIA GPU 数量
        info["cores"] = pynvml.nvmlDeviceGetCount()
        # 遍历每个 GPU，获取 GPU 信息
        for i in range(info["cores"]):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            # 获取 GPU 型号
            info["type"].append(pynvml.nvmlDeviceGetName(handle))

    except pynvml.NVMLError as e:
        print(f"Error: {e}")
    finally:
        # 结束 NVML
        pynvml.nvmlShutdown()
        return info


def __get_command():
    """获取执行训练时的完整命令行信息
    比如在运行`python main.py -i 123`时，full_command为`main.py -i 123`
    """
    import sys

    full_command = " ".join(sys.argv)
    return full_command


def __get_pip_requirement():
    """获取当前项目下的全部Python环境，是1个很长的、带换行的文件列表，建议后续存储在swanlog目录下"""
    try:
        # 运行pip命令获取当前环境下的环境目录
        result = subprocess.run(["pip", "freeze"], stdout=subprocess.PIPE, text=True)

        # 检查命令是否成功运行
        if result.returncode == 0:
            return result.stdout
        else:
            print(f"Error: {result.stderr}")
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def __get_pip_requirement():
    """获取当前项目下的全部Python环境，是1个很长的、带换行的文件列表，建议后续存储在swanlog目录下"""
    try:
        # 运行pip命令获取当前环境下的环境目录
        result = subprocess.run(["pip", "freeze"], stdout=subprocess.PIPE, text=True)

        # 检查命令是否成功运行
        if result.returncode == 0:
            return result.stdout
        else:
            print(f"Error: {result.stderr}")
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def get_system_info():
    """获取系统信息"""
    return {
        "hostname": socket.gethostname(),
        "os": platform.platform(),
        "python": platform.python_version(),
        "executable": sys.executable,  # python 解释器路径
        "git_remote": __get_remote_url(),  # 获取远程仓库的链接
        "cpu": multiprocessing.cpu_count(),  # cpu 核心数
        "gpu": __get_gpu_info(),  # gpu 相关信息
    }
