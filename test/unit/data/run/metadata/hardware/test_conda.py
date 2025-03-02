"""
@author: cunyue
@file: test_conda.py
@time: 2025/3/2 13:08
@description: 测试conda
"""

import subprocess

import pytest
import yaml

from swanlab.data.run.metadata.conda import get_conda

has_conda = subprocess.run(["conda", "--version"], shell=True, capture_output=True).returncode == 0


@pytest.mark.skipif(not has_conda, reason="conda is not installed")
def test_conda():
    conda_info = get_conda()
    assert isinstance(conda_info, str)
    # 可被yaml解析
    load_conda = yaml.safe_load(conda_info)
    assert load_conda is not None
    assert isinstance(load_conda, dict)
