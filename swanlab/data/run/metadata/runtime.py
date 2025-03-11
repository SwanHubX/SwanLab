"""
@author: cunyue
@file: runtime.py
@time: 2024/11/18 15:13
@description: 运行时信息，操作系统信息等
"""

import os
import platform
import socket
import subprocess
import sys


def get_runtime_info():
    return {
        **get_computer_info(),
        **get_python_info(),
        **get_git_info(),
    }


def get_computer_info():
    return {
        "os": platform.platform(),
        "hostname": socket.gethostname(),
        "pid": os.getpid(),
        "cwd": os.getcwd(),
    }


def get_python_info():
    return {
        "python": platform.python_version(),
        "python_verbose": sys.version,
        "executable": sys.executable,
        "command": " ".join(sys.argv),
    }


# ---------------------------------- git信息 ----------------------------------


def get_git_info():
    """获取git信息"""
    return {
        "git_remote": get_remote_url(),
        "git_info": get_git_branch_and_commit(),
    }


def get_remote_url():
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
            url = result.stdout.strip()
            return parse_git_url(url)
        else:
            return None
    except Exception as e:  # noqa
        return None

def parse_git_url(url):
    """Return the remote URL of a git repository."""
    if url.startswith("git@"):
        parts = url[4:].split("/", 1)
        host, path = parts[0], parts[1] if len(parts) > 1 else ""
        if ":" in host:
            host, port = host.rsplit(":", 1)
            url = f"https://{host}:{port}/{path}" if port.isdigit() else f"https://{host}/{port}/{path}"
        else:
            url = f"https://{host}/{path}"
    return url[:-4] if url.endswith(".git") else url

def replace_second_colon(input_string, replacement):
    """Replace the second colon in a string."""
    first_colon = input_string.find(":")
    second_colon = input_string.find(":", first_colon + 1) if first_colon != -1 else -1
    return input_string[:second_colon] + replacement + input_string[second_colon + 1:] if second_colon != -1 else input_string


def get_git_branch_and_commit():
    """获取项目git的分支名和该分支下最新提交的hash"""
    try:
        # 获取当前分支名称
        branch_process = subprocess.Popen(
            ["git", "branch", "--show-current"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        branch_name, branch_err = branch_process.communicate()
        if branch_process.returncode != 0:
            return None, None

        branch_name = branch_name.strip()

        # 获取当前分支的最新提交hash
        commit_process = subprocess.Popen(
            ["git", "rev-parse", branch_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        commit_hash, commit_err = commit_process.communicate()

        # 如果无法获取最新提交hash, 那么就返回branch_name和None
        if commit_process.returncode != 0:
            return branch_name, None

        commit_hash = commit_hash.strip()
        return branch_name, commit_hash

    except Exception as e:  # noqa
        return None, None
