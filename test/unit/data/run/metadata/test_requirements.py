"""
@author: cunyue
@file: test_requirements.py
@time: 2025/4/25 22:09
@description: 测试获取python依赖
"""

from swanlab.data.run.metadata.requirements import get_requirements


def test_get_requirements():
    """
    测试获取python依赖
    正常情况下有python都应该获取成功
    """
    requirements = get_requirements()
    assert requirements is not None
    assert isinstance(requirements, str)
