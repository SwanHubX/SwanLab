"""
@author: cunyue
@file: test_conda.py
@time: 2025/3/2 13:08
@description: 测试conda
"""

import subprocess

from swanlab.data.run.metadata.conda import get_conda


has_conda = subprocess.run(["conda", "--version"], shell=True, capture_output=True).returncode == 0

def test_conda():
    conda_info = get_conda()
    if not has_conda:
        assert conda_info is None
        return
    assert isinstance(conda_info, dict)
