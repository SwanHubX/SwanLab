"""
@author: cunyue
@file: conftest.py
@time: 2026/3/12
@description: benchmark 公共 fixtures 与收集约束

``e2e/bench_e2e.py`` 与 ``multi_process/bench_multi_process.py`` 是纯命令行驱动（没有
``test_*`` 入口，单跑就是百万级场景），不参与 pytest 收集；其余 ``bench_*.py`` 保持原有
``make bench`` 行为不变。
"""

# 只排除无测试入口的脚本与产物目录，保留历史 bench_*.py 的 pytest 入口。
collect_ignore_glob = [
    "e2e/bench_e2e.py",
    "multi_process/*",
    "mock_cloud/*",
    "results/*",
]
