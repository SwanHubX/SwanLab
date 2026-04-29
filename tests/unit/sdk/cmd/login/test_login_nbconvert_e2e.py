"""
@author: cunyue
@file: test_login_nbconvert_e2e.py
@time: 2026/4/29
@description: 测试 swanlab.login() 在 jupyter nbconvert --execute 环境下的行为 (issue #1334)

在 nbconvert --execute 环境中，kernel 不支持 stdin（_allow_stdin=False），
调用 swanlab.login() 不带参数时应抛出 AuthenticationError 而非 StdinNotImplementedError。
"""

import json
import os
import subprocess
import sys

import pytest

NOTEBOOK_TEMPLATE = {
    "cells": [
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import os\n",
                "os.environ.pop('SWANLAB_API_KEY', None)\n",
                "from swanlab.exceptions import AuthenticationError\n",
                "from swanlab.sdk.cmd.login import login\n",
                "try:\n",
                "    login()\n",
                "    print('LOGIN_SUCCESS')\n",
                "except AuthenticationError:\n",
                "    print('AuthenticationError_CAUGHT')\n",
                "except Exception as e:\n",
                "    print(f'UNEXPECTED_{type(e).__name__}')\n",
            ],
        }
    ],
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.11.0"},
    },
    "nbformat": 4,
    "nbformat_minor": 4,
}


def _execute_nbconvert(notebook_path: str, home_dir: str) -> tuple:
    """在隔离 HOME 下执行 nbconvert --execute，返回 (subprocess result, cell stdout)"""
    output_path = notebook_path.replace(".ipynb", "_executed.ipynb")
    env = os.environ.copy()
    env["HOME"] = home_dir
    env.pop("SWANLAB_API_KEY", None)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "jupyter",
            "nbconvert",
            "--execute",
            notebook_path,
            "--to",
            "notebook",
            "--output",
            output_path,
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    # 读取执行后的 notebook
    with open(output_path) as f:
        nb = json.load(f)
    # 提取第一个 code cell 的 stdout
    stdout_text = ""
    for out in nb["cells"][0].get("outputs", []):
        if out.get("output_type") == "stream" and out.get("name") == "stdout":
            stdout_text = "".join(out.get("text", []))
    return result, stdout_text


@pytest.mark.skipif(
    not subprocess.run([sys.executable, "-m", "jupyter", "--version"], capture_output=True).returncode == 0,
    reason="jupyter not installed",
)
def test_login_raises_authentication_error_in_nbconvert(tmp_path):
    """在 nbconvert --execute 环境下，login() 无 API key 应抛出 AuthenticationError"""
    # 写入 notebook
    nb_path = str(tmp_path / "test_login.ipynb")
    with open(nb_path, "w") as f:
        json.dump(NOTEBOOK_TEMPLATE, f)
    # 隔离 HOME，确保无 .netrc 凭证
    home_dir = str(tmp_path / "home")
    os.makedirs(home_dir, exist_ok=True)

    result, stdout = _execute_nbconvert(nb_path, home_dir)
    assert "AuthenticationError_CAUGHT" in stdout, (
        f"Expected AuthenticationError_CAUGHT, got: {stdout}\nstderr: {result.stderr}"
    )
