"""
@author: cunyue
@file: test_helper.py
@time: 2026/3/14
@description: 测试 swanlab.sdk.cmd.helper 中的 with_cmd_lock 装饰器
"""

import threading
import time

from swanlab.sdk.cmd.helper import with_cmd_lock


class TestWithCmdLock:
    def test_wraps_preserves_metadata(self):
        """@wraps 应保留被装饰函数的 __name__ 和 __doc__"""

        def my_func():
            """My docstring"""

        wrapped = with_cmd_lock(my_func)

        assert wrapped.__name__ == "my_func"
        assert wrapped.__doc__ == "My docstring"

    def test_lock_serializes_concurrent_calls(self):
        """并发调用两个装饰函数时，执行顺序不得交错（严格串行化）"""
        events = []

        @with_cmd_lock
        def record_event(name):
            events.append(f"{name}:start")
            time.sleep(0.03)
            events.append(f"{name}:end")

        t1 = threading.Thread(target=record_event, args=("A",))
        t2 = threading.Thread(target=record_event, args=("B",))

        t1.start()
        time.sleep(0.005)  # 确保 t1 先抢到锁
        t2.start()

        t1.join()
        t2.join()

        # 无论谁先执行，各自的 start/end 必须紧邻，不允许交错
        assert events in [
            ["A:start", "A:end", "B:start", "B:end"],
            ["B:start", "B:end", "A:start", "A:end"],
        ]
