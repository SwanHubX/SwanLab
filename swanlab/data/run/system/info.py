#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024-03-19 15:37:02
@File: swanlab/data/system/info.py
@IDE: vscode
@Description:
    获取系统信息
"""
import platform
import socket
import sys
import subprocess
import multiprocessing
import pynvml
from swanlab.log import swanlog
import os


def __replace_second_colon(input_string, replacement):
    """替换字符串中第二个‘:’"""
    first_colon_index = input_string.find(":")
    if first_colon_index != -1:
        second_colon_index = input_string.find(":", first_colon_index + 1)
        if second_colon_index != -1:
            return input_string[:second_colon_index] + replacement + input_string[second_colon_index + 1:]
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
            swanlog.debug(f"An error occurred: {result.stderr}")
            return None
    except Exception as e:
        swanlog.debug(f"An error occurred: {e}")
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
            swanlog.debug("Error getting branch name:", branch_err)
            return None, None

        branch_name = branch_name.strip()

        # 获取当前分支的最新提交hash
        commit_process = subprocess.Popen(
            ["git", "rev-parse", branch_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        commit_hash, commit_err = commit_process.communicate()

        # 如果无法获取最新提交hash, 那么就返回branch_name和None
        if commit_process.returncode != 0:
            swanlog.debug("Error getting commit hash:", commit_err)
            return branch_name, None

        commit_hash = commit_hash.strip()
        return branch_name, commit_hash

    except Exception as e:
        swanlog.debug(f"An error occurred: {e}")
        return None, None


def __get_nvidia_gpu_info():
    """获取 GPU 信息"""
    info = {"cores": None, "type": [], "memory": []}
    try:
        pynvml.nvmlInit()
    except:
        return None

    try:
        # 获取 NVIDIA GPU 数量
        info["cores"] = pynvml.nvmlDeviceGetCount()
        # 遍历每个 GPU，获取 GPU 信息
        for i in range(info["cores"]):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            # 获取 GPU 型号
            gpu_name = pynvml.nvmlDeviceGetName(handle)  # types: bytes | str
            if isinstance(gpu_name, bytes):  # Fix for pynvml 早期版本，关联 issue: #605
                gpu_name = gpu_name.decode("utf-8")
            info["type"].append(gpu_name)
            # 获取 GPU 的总显存, 单位为GB
            info["memory"].append(round(pynvml.nvmlDeviceGetMemoryInfo(handle).total / (1024 ** 3)))

    except pynvml.NVMLError as e:
        swanlog.debug(f"An error occurred when getting GPU info: {e}")
        pass
    finally:
        # 结束 NVML
        pynvml.nvmlShutdown()
        return info


def __get_apple_gpu_info():
    import json

    info = {"cores": None, "type": [], "memory": []}

    # 使用system_profiler命令以JSON格式获取GPU信息
    try:
        result = subprocess.run(["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True)
        gpu_name = json.loads(result.stdout)["SPHardwareDataType"][0]["chip_type"]
        memory = json.loads(result.stdout)["SPHardwareDataType"][0]["physical_memory"]
        memory = str(memory).lower().replace("gb", "")
        number_processors = json.loads(result.stdout)["SPHardwareDataType"][0]["number_processors"]
    except:
        return None

    info["type"].append(gpu_name)
    info["memory"].append(memory)
    info["cores"] = number_processors

    # TODO: Apple设备硬件监控时再解除注释，二进制文件apple_gpu_stats来自https://github.com/wandb/wandb/blob/main/wandb/bin/apple_gpu_stats
    #
    # MAX_POWER_WATTS = 16.5
    #
    # import pathlib
    # binary_path = (pathlib.Path(sys.modules["swanlab"].__path__[0]) / "bin" / "apple_gpu_stats").resolve()
    # try:
    #     command = [str(binary_path), "--json"]
    #     output = (subprocess.check_output(command, universal_newlines=True).strip().split("\n"))[0]
    #     raw_stats = json.loads(output)
    #     stats = {
    #         "gpu": raw_stats["utilization"],
    #         "memoryAllocated": raw_stats["mem_used"],
    #         "temp": raw_stats["temperature"],
    #         "powerWatts": raw_stats["power"],
    #         "powerPercent": (raw_stats["power"] / MAX_POWER_WATTS) * 100,
    #     }
    # except:
    #     swanlog.debug(f"Apple GPU stats failed to obtain.")

    return info


def __get_gpu_info():
    gpu_info = __get_nvidia_gpu_info()
    if gpu_info is not None:
        return gpu_info

    apple_info = __get_apple_gpu_info()
    if apple_info is not None:
        return apple_info

    return {"cores": None, "type": [], "memory": []}


def __get_command():
    """获取执行训练时的完整命令行信息
    比如在运行`python main.py -i 123`时，full_command为`main.py -i 123`
    """
    import sys

    full_command = " ".join(sys.argv)
    return full_command


def __get_memory_size():
    """获取内存大小"""
    import psutil

    try:
        # 获取系统总内存大小
        mem = psutil.virtual_memory()
        total_memory = round(mem.total / (1024 ** 3))  # 单位为GB
        return total_memory
    except Exception as e:
        swanlog.debug(f"An error occurred when getting memory size: {e}")
        return None


def __get_cwd():
    """获取当前工作目录路径"""
    import os

    try:
        cwd = os.getcwd()
        return cwd
    except Exception as e:
        swanlog.debug(f"An error occurred when getting current working directory: {e}")
        return None


def get_requirements() -> str:
    """获取当前项目下的全部Python环境，是1个很长的、带换行的文件列表，建议后续存储在swanlog目录下"""
    try:
        # 运行pip命令获取当前环境下的环境目录
        result = subprocess.run(["pip", "list", "--format=freeze"], stdout=subprocess.PIPE, text=True)

        # 检查命令是否成功运行
        if result.returncode == 0:
            return result.stdout
        else:
            swanlog.debug(f"An error occurred when getting requirements:{result.stderr}")
            return None
    except Exception as e:
        swanlog.debug(f"An error occurred: {e}")
        return None


def get_system_info(version: str, logdir: str):
    """获取系统信息
    :param version: swanlab版本号
    :param logdir: swanlab日志目录
    """
    return {
        "swanlab": {"version": version, "logdir": logdir},  # swanlab 版本号和日志目录
        "hostname": socket.gethostname(),  # 主机名
        "os": platform.platform(),  # 操作系统
        "python": platform.python_version(),  # python版本
        "python_verbose": sys.version,  # python详细版本
        "executable": sys.executable,  # python 解释器路径
        "git_remote": __get_remote_url(),  # 获取远程仓库的链接
        "cpu": multiprocessing.cpu_count(),  # cpu 核心数
        "gpu": __get_gpu_info(),  # gpu 相关信息
        "git_info": __get_git_branch_and_commit(),  # git 分支和最新 commit 信息
        "command": __get_command(),  # 完整命令行信息
        "memory": __get_memory_size(),  # 内存大小
        "cwd": __get_cwd(),  # 当前工作目录路径
        "pid": os.getpid(),  # 当前进程ID
    }
